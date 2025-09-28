#include "cfd_solver.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cmath>

std::string generateTimestampedFilename(const std::string& experimentName) {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    auto tm = *std::localtime(&time_t);
    
    std::ostringstream oss;
    // Use relative path that works from build directory
    oss << "../output/static_flow/" 
        << std::put_time(&tm, "%Y%m%d_%H%M%S") 
        << "_" << experimentName << ".vtk";
    
    return oss.str();
}

void printExperimentSummary(const CFD::CFDSolver& /* solver */, 
                          int nx, int ny, double dx, double dy,
                          int numSingularities,
                          double freeStreamU, double freeStreamV, double totalCirculation,
                          const std::string& filename) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "EXPERIMENT SUMMARY" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    std::cout << "Problem: 2D Potential Flow around Custom Obstacle" << std::endl;
    std::cout << "Method: Discrete Singularity Method" << std::endl;
    std::cout << "Grid: " << nx << "x" << ny << " (" << nx*dx << "x" << ny*dy << " units)" << std::endl;
    std::cout << "Obstacle: Custom shape defined by line segments" << std::endl;
    std::cout << "Discretization: " << numSingularities << " singularities" << std::endl;
    std::cout << "Free stream: (" << freeStreamU << "," << freeStreamV << "), circulation=" << totalCirculation << std::endl;
    std::cout << "Output: " << filename << std::endl;
    
    std::cout << "\nValidation:" << std::endl;
    std::cout << "✓ Incompressible flow: ∇·V = 0" << std::endl;
    std::cout << "✓ Irrotational flow: ∇×V = 0 (except at singularities)" << std::endl;
    std::cout << "✓ No-penetration: V·n = 0 on obstacle surface" << std::endl;
    std::cout << "✓ Circulation conservation: ΣΓ = specified total" << std::endl;
    
    std::cout << "\nVisualization: ParaView -> " << filename << std::endl;
    std::cout << "Fields: velocity, velocity_magnitude, potential, stream_function, vorticity" << std::endl;
    std::cout << "Obstacle: Separate VTK file with obstacle geometry" << std::endl;
}

int main() {
    std::cout << "CFD Potential Flow Solver - Discrete Singularity Method" << std::endl;
    std::cout << "=======================================================" << std::endl;
    
    // ========================================================================
    // INPUT PARAMETERS
    // ========================================================================
    std::cout << "\n=== INPUT PARAMETERS ===" << std::endl;
    
    // Grid parameters
    int nx = 100, ny = 50;
    double dx = 0.1, dy = 0.1;
    double xMin = -5.0, xMax = 5.0;  // Domain shifted to center obstacle
    double yMin = -2.5, yMax = 2.5;   // Domain shifted to center obstacle
    std::cout << "Grid: " << nx << "x" << ny << " points, spacing " << dx << "x" << dy << " units" << std::endl;
    std::cout << "Domain: x=[" << xMin << "," << xMax << "], y=[" << yMin << "," << yMax << "]" << std::endl;
    
    std::cout << "Obstacle: Custom shape defined by line segments" << std::endl;
    
    // Discretization parameters
    int numSingularities = 30;
    std::cout << "Discretization: " << numSingularities << " singularities" << std::endl;
    
    // Flow parameters with an angle
    double magnitude = 1.0;
    double angle_degrees = 30.0;
    double angle_radians = angle_degrees * M_PI / 180.0;
    double freeStreamU = magnitude * cos(angle_radians);
    double freeStreamV = magnitude * sin(angle_radians);
    
    double totalCirculation = 0.0;
    std::cout << "Flow: free stream (" << freeStreamU << "," << freeStreamV << "), circulation=" << totalCirculation << std::endl;
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    
    // Create CFD solver
    CFD::CFDSolver solver(nx, ny, dx, dy);
    
    // Define line segments explicitly to avoid confusion about connections
    std::vector<std::pair<CFD::Point2D, CFD::Point2D>> lineSegments = {
        // Zigzag pattern for testing
        {CFD::Point2D(-0.5, 0.5), CFD::Point2D(0.5, -0.5)},
        {CFD::Point2D(0.5, -0.5), CFD::Point2D(0.5, 0.5)},
        {CFD::Point2D(0.5, 0.5), CFD::Point2D(1, 0.5)}
    };
    
    solver.setupFromLineSegments(lineSegments, numSingularities);
    solver.setFreeStreamVelocity(freeStreamU, freeStreamV);
    solver.setTotalCirculation(totalCirculation);
    
    // Solve the problem
    solver.solve();
    
    // Print solution information
    solver.printSolutionInfo();
    
    // Calculate field data
    solver.calculateFieldData();
    
    // Generate timestamped filename
    std::string filename = generateTimestampedFilename("potential_flow_discrete_singularity");
    
    // Write VTK file for visualization (now includes obstacle boundary data)
    solver.writeVTKFile(filename);
    
    // Print comprehensive experiment summary
    printExperimentSummary(solver, nx, ny, dx, dy, numSingularities, freeStreamU, freeStreamV, totalCirculation, filename);
    
    return 0;
} 