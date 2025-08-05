#include "cfd_solver.h"
#include "complex_math.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <sstream>

namespace CFD {

// FieldData implementation
FieldData::FieldData(int nx, int ny, double dx, double dy) 
    : nx(nx), ny(ny), dx(dx), dy(dy) {
    
    u.resize(nx, std::vector<double>(ny));
    v.resize(nx, std::vector<double>(ny));
    potential.resize(nx, std::vector<double>(ny));
    streamFunction.resize(nx, std::vector<double>(ny));
    velocityMagnitude.resize(nx, std::vector<double>(ny));
    vorticity.resize(nx, std::vector<double>(ny));
}

// CFDSolver implementation
CFDSolver::CFDSolver(int nx, int ny, double dx, double dy) 
    : freeStreamVelocity(1.0, 0.0), totalCirculation(0.0), fieldData(nx, ny, dx, dy) {
    
    std::cout << "CFD Solver: " << nx << "x" << ny << " grid, domain " 
              << nx*dx << "x" << ny*dy << " units" << std::endl;
    
    // Validate grid parameters
    if (nx <= 0 || ny <= 0) {
        std::cerr << "ERROR: Grid dimensions must be positive!" << std::endl;
        throw std::invalid_argument("Invalid grid dimensions");
    }
    if (dx <= 0 || dy <= 0) {
        std::cerr << "ERROR: Grid spacing must be positive!" << std::endl;
        throw std::invalid_argument("Invalid grid spacing");
    }
}



void CFDSolver::setupFromLineSegments(const std::vector<std::pair<Point2D, Point2D>>& lineSegments, int numSingularities) {
    std::cout << "Obstacle: custom shape with " << lineSegments.size() << " line segments, " 
              << numSingularities << " singularities" << std::endl;
    
    // Validate obstacle parameters
    if (lineSegments.empty()) {
        std::cerr << "ERROR: Need at least 1 line segment!" << std::endl;
        throw std::invalid_argument("Invalid line segments");
    }
    if (numSingularities < 2) {
        std::cerr << "ERROR: Need at least 2 singularities!" << std::endl;
        throw std::invalid_argument("Invalid number of singularities");
    }
    
    obstacle = ObstacleContour::createFromLineSegments(lineSegments, numSingularities);
    obstacle.generateDiscretization(numSingularities);
    singularities = obstacle.getSingularities();
    collocationPoints = obstacle.getCollocationPoints();
}

void CFDSolver::setFreeStreamVelocity(double u, double v) {
    double magnitude = std::sqrt(u*u + v*v);
    double angle = std::atan2(v, u) * 180.0 / M_PI;
    std::cout << "Free stream: (" << u << "," << v << ") = " << magnitude 
              << " units at " << angle << "°" << std::endl;
    freeStreamVelocity = Point2D(u, v);
}

void CFDSolver::setTotalCirculation(double circulation) {
    std::cout << "Total circulation: " << circulation << std::endl;
    totalCirculation = circulation;
}

void CFDSolver::solve() {
    std::cout << "\n=== SOLVING LINEAR SYSTEM ===" << std::endl;
    std::cout << "System: " << singularities.size() << "x" << singularities.size() 
              << " (" << collocationPoints.size() << " collocation + 1 circulation)" << std::endl;
    
    // Validate system setup
    if (singularities.size() != collocationPoints.size() + 1) {
        std::cerr << "ERROR: System is underdetermined!" << std::endl;
        throw std::runtime_error("Invalid system setup");
    }
    
    // Create and solve the linear system
    LinearSystemSolver solver(singularities.size());
    solver.buildSystem(singularities, collocationPoints, freeStreamVelocity, totalCirculation);
    
    // Solve for circulations
    std::vector<double> circulations = solver.solve();
    
    // Update singularities with solved circulations
    for (size_t i = 0; i < singularities.size(); ++i) {
        singularities[i].circulation = circulations[i];
    }
    
    // Log solution statistics
    double minCirculation = *std::min_element(circulations.begin(), circulations.end());
    double maxCirculation = *std::max_element(circulations.begin(), circulations.end());
    double sumCirculations = std::accumulate(circulations.begin(), circulations.end(), 0.0);
    
    std::cout << "Circulation range: [" << std::fixed << std::setprecision(3) 
              << minCirculation << ", " << maxCirculation << "], sum=" << sumCirculations << std::endl;
}

void CFDSolver::calculateFieldData() {
    std::cout << "\n=== CALCULATING FIELD DATA ===" << std::endl;
    std::cout << "Computing " << fieldData.nx * fieldData.ny << " grid points..." << std::endl;
    
    // Statistics for field data
    double minVelocity = std::numeric_limits<double>::max();
    double maxVelocity = std::numeric_limits<double>::lowest();
    double minPotential = std::numeric_limits<double>::max();
    double maxPotential = std::numeric_limits<double>::lowest();
    double minStreamFunction = std::numeric_limits<double>::max();
    double maxStreamFunction = std::numeric_limits<double>::lowest();
    double minVorticity = std::numeric_limits<double>::max();
    double maxVorticity = std::numeric_limits<double>::lowest();
    
    int pointsNearSingularities = 0;
    
    // Calculate domain bounds for shifted coordinate system
    double xMin = -5.0;
    double yMin = -2.5;
    
    for (int i = 0; i < fieldData.nx; ++i) {
        for (int j = 0; j < fieldData.ny; ++j) {
            double x = xMin + i * fieldData.dx;
            double y = yMin + j * fieldData.dy;
            
            // Initialize with free stream
            double u = freeStreamVelocity.x;
            double v = freeStreamVelocity.y;
            double potential = x * freeStreamVelocity.x + y * freeStreamVelocity.y;
            double streamFunction = y * freeStreamVelocity.x - x * freeStreamVelocity.y;
            
            // Add contributions from all singularities
            for (const auto& singularity : singularities) {
                // Check if point is very close to singularity
                double distToSingularity = std::sqrt((x - singularity.position.x)*(x - singularity.position.x) + 
                                                   (y - singularity.position.y)*(y - singularity.position.y));
                if (distToSingularity < 0.1) {
                    pointsNearSingularities++;
                }
                
                // Velocity contribution
                Point2D vortexVel = vortexVelocity(x, y, 
                    singularity.position.x, singularity.position.y, 
                    singularity.circulation);
                u += vortexVel.x;
                v += vortexVel.y;
                
                // Potential contribution (equation 7.1.22')
                double theta = angularFunction(x, y, 
                    singularity.position.x, singularity.position.y);
                potential += singularity.circulation * theta;
                
                // Stream function contribution (equation 7.1.23')
                double psi = vortexStreamFunction(x, y, 
                    singularity.position.x, singularity.position.y, 
                    singularity.circulation);
                streamFunction += psi;
            }
            
            // Store field data
            fieldData.u[i][j] = u;
            fieldData.v[i][j] = v;
            fieldData.potential[i][j] = potential;
            fieldData.streamFunction[i][j] = streamFunction;
            fieldData.velocityMagnitude[i][j] = std::sqrt(u*u + v*v);
            
            // Update statistics
            minVelocity = std::min(minVelocity, fieldData.velocityMagnitude[i][j]);
            maxVelocity = std::max(maxVelocity, fieldData.velocityMagnitude[i][j]);
            minPotential = std::min(minPotential, potential);
            maxPotential = std::max(maxPotential, potential);
            minStreamFunction = std::min(minStreamFunction, streamFunction);
            maxStreamFunction = std::max(maxStreamFunction, streamFunction);
            
            // Calculate vorticity (simplified)
            if (i > 0 && i < fieldData.nx-1 && j > 0 && j < fieldData.ny-1) {
                double dudy = (fieldData.u[i][j+1] - fieldData.u[i][j-1]) / (2.0 * fieldData.dy);
                double dvdx = (fieldData.v[i+1][j] - fieldData.v[i-1][j]) / (2.0 * fieldData.dx);
                fieldData.vorticity[i][j] = dvdx - dudy;
                minVorticity = std::min(minVorticity, fieldData.vorticity[i][j]);
                maxVorticity = std::max(maxVorticity, fieldData.vorticity[i][j]);
            } else {
                fieldData.vorticity[i][j] = 0.0;
            }
        }
    }
    
    std::cout << "Field statistics:" << std::endl;
    std::cout << "  Velocity: [" << std::fixed << std::setprecision(3) 
              << minVelocity << ", " << maxVelocity << "]" << std::endl;
    std::cout << "  Potential: [" << std::fixed << std::setprecision(3) 
              << minPotential << ", " << maxPotential << "]" << std::endl;
    std::cout << "  Stream function: [" << std::fixed << std::setprecision(3) 
              << minStreamFunction << ", " << maxStreamFunction << "]" << std::endl;
    std::cout << "  Vorticity: [" << std::fixed << std::setprecision(3) 
              << minVorticity << ", " << maxVorticity << "]" << std::endl;
    std::cout << "  Points near singularities: " << pointsNearSingularities << std::endl;
}

void CFDSolver::writeVTKFile(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    // Write VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "CFD 2D Potential Flow - Discrete Singularity Method\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << fieldData.nx << " " << fieldData.ny << " 1\n";
    file << "POINTS " << fieldData.nx * fieldData.ny << " float\n";
    
    // Write grid points with shifted coordinates
    double xMin = -5.0, yMin = -2.5;
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            double x = xMin + i * fieldData.dx;
            double y = yMin + j * fieldData.dy;
            file << x << " " << y << " 0.0\n";
        }
    }
    
    // Write velocity field
    file << "\nPOINT_DATA " << fieldData.nx * fieldData.ny << "\n";
    file << "VECTORS velocity float\n";
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            file << fieldData.u[i][j] << " " << fieldData.v[i][j] << " 0.0\n";
        }
    }
    
    // Write velocity magnitude
    file << "\nSCALARS velocity_magnitude float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            file << fieldData.velocityMagnitude[i][j] << "\n";
        }
    }
    
    // Write potential
    file << "\nSCALARS potential float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            file << fieldData.potential[i][j] << "\n";
        }
    }
    
    // Write stream function
    file << "\nSCALARS stream_function float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            file << fieldData.streamFunction[i][j] << "\n";
        }
    }
    
    // Write vorticity
    file << "\nSCALARS vorticity float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < fieldData.ny; ++j) {
        for (int i = 0; i < fieldData.nx; ++i) {
            file << fieldData.vorticity[i][j] << "\n";
        }
    }
    
    file.close();
    std::cout << "VTK file written: " << filename << std::endl;
}



void CFDSolver::writeThickObstacleVTKFile(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    // Get obstacle vertices from the obstacle contour
    const auto& vertices = obstacle.getVertices();
    
    if (vertices.empty()) {
        std::cerr << "Error: No obstacle vertices available!" << std::endl;
        return;
    }
    
    // Create a thick line representation by creating rectangles along each segment
    double thickness = 0.025;  // Thickness of the obstacle line (half the previous width)
    int numSegments = vertices.size();
    
    // Calculate total number of points (4 points per segment for rectangle)
    int totalPoints = numSegments * 4;
    
    // Write VTK header for unstructured grid
    file << "# vtk DataFile Version 3.0\n";
    file << "CFD 2D Thick Obstacle Geometry\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write points (4 points per segment to form rectangles)
    file << "POINTS " << totalPoints << " float\n";
    
    // Create segments based on the actual obstacle geometry
    // Since we now use line segments, we can create thick rectangles for each segment
    std::vector<std::pair<size_t, size_t>> segments;
    
    // For now, create segments between consecutive vertices
    // This will work for the line segments approach
    for (size_t i = 0; i < vertices.size() - 1; ++i) {
        segments.push_back({i, i + 1});
    }
    
    for (const auto& segment : segments) {
        size_t i = segment.first;
        size_t next = segment.second;
        
        // Calculate segment direction
        double dx = vertices[next].x - vertices[i].x;
        double dy = vertices[next].y - vertices[i].y;
        double length = std::sqrt(dx*dx + dy*dy);
        
        if (length < 1e-10) continue;  // Skip zero-length segments
        
        // Normalize direction
        dx /= length;
        dy /= length;
        
        // Calculate perpendicular direction (normal)
        double nx = -dy;
        double ny = dx;
        
        // Create rectangle points for this segment
        double x1 = vertices[i].x + nx * thickness;
        double y1 = vertices[i].y + ny * thickness;
        double x2 = vertices[i].x - nx * thickness;
        double y2 = vertices[i].y - ny * thickness;
        double x3 = vertices[next].x - nx * thickness;
        double y3 = vertices[next].y - ny * thickness;
        double x4 = vertices[next].x + nx * thickness;
        double y4 = vertices[next].y + ny * thickness;
        
        // Write the 4 points of this rectangle
        file << x1 << " " << y1 << " 0.0\n";
        file << x2 << " " << y2 << " 0.0\n";
        file << x3 << " " << y3 << " 0.0\n";
        file << x4 << " " << y4 << " 0.0\n";
    }
    
    // Write cells (one quad per segment)
    int totalCells = segments.size();
    int totalCellListSize = totalCells * 5;  // 4 points + 1 count per cell
    
    file << "\nCELLS " << totalCells << " " << totalCellListSize << "\n";
    for (int i = 0; i < totalCells; ++i) {
        file << "4";  // Number of points in this cell (quad)
        for (int j = 0; j < 4; ++j) {
            file << " " << (i * 4 + j);  // Point indices for this quad
        }
        file << "\n";
    }
    
    // Write cell types (all are quads)
    file << "\nCELL_TYPES " << totalCells << "\n";
    for (int i = 0; i < totalCells; ++i) {
        file << "9\n";  // VTK_QUAD
    }
    
    // Add cell data for visualization
    file << "\nCELL_DATA " << totalCells << "\n";
    file << "SCALARS obstacle_id float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < totalCells; ++i) {
        file << "1.0\n";  // All cells belong to the obstacle
    }
    
    file.close();
    std::cout << "Thick obstacle VTK file written: " << filename << std::endl;
}

void CFDSolver::printSolutionInfo() const {
    std::cout << "\n=== Solution Information ===" << std::endl;
    std::cout << "Free stream velocity: (" << freeStreamVelocity.x << ", " << freeStreamVelocity.y << ")" << std::endl;
    std::cout << "Total circulation: " << totalCirculation << std::endl;
    std::cout << "Number of singularities: " << singularities.size() << std::endl;
    
    std::cout << "\nSingularity circulations:" << std::endl;
    for (size_t i = 0; i < singularities.size(); ++i) {
        const auto& s = singularities[i];
        std::cout << "  Γ[" << i << "] = " << s.circulation 
                  << " at (" << s.position.x << ", " << s.position.y << ")" << std::endl;
    }
    
    double total = 0.0;
    for (const auto& s : singularities) {
        total += s.circulation;
    }
    std::cout << "Sum of circulations: " << total << std::endl;
}

} // namespace CFD 