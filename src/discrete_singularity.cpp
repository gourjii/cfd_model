#include "discrete_singularity.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric> // Required for std::accumulate

namespace CFD {



ObstacleContour ObstacleContour::createFromLineSegments(const std::vector<std::pair<Point2D, Point2D>>& lineSegments, int numPoints) {
    ObstacleContour contour;
    
    if (lineSegments.empty()) {
        std::cerr << "ERROR: Need at least 1 line segment!" << std::endl;
        return contour;
    }
    
    std::cout << "Creating obstacle from " << lineSegments.size() << " line segments:" << std::endl;
    for (size_t i = 0; i < lineSegments.size(); ++i) {
        const auto& segment = lineSegments[i];
        std::cout << "  Segment " << i << ": (" << segment.first.x << ", " << segment.first.y 
                  << ") to (" << segment.second.x << ", " << segment.second.y << ")" << std::endl;
    }
    
    // Calculate total perimeter length
    double totalLength = 0.0;
    for (const auto& segment : lineSegments) {
        double dx = segment.second.x - segment.first.x;
        double dy = segment.second.y - segment.first.y;
        totalLength += std::sqrt(dx*dx + dy*dy);
    }
    
    // Distribute points along each segment based on its length
    for (const auto& segment : lineSegments) {
        // Calculate segment length
        double dx = segment.second.x - segment.first.x;
        double dy = segment.second.y - segment.first.y;
        double segmentLength = std::sqrt(dx*dx + dy*dy);
        
        // Calculate number of points for this segment (proportional to length)
        int segmentPoints = std::max(2, (int)(numPoints * segmentLength / totalLength));
        
        // Add points along this segment (excluding the last point to avoid duplicates)
        for (int j = 0; j < segmentPoints - 1; ++j) {
            double t = j / (segmentPoints - 1.0);
            double x = segment.first.x + (segment.second.x - segment.first.x) * t;
            double y = segment.first.y + (segment.second.y - segment.first.y) * t;
            contour.vertices.push_back(Point2D(x, y));
        }
        
        // Add the end point of this segment
        contour.vertices.push_back(segment.second);
    }
    
    // Remove duplicate vertices (keep only unique ones)
    if (!contour.vertices.empty()) {
        std::vector<Point2D> uniqueVertices;
        uniqueVertices.push_back(contour.vertices[0]);
        
        for (size_t i = 1; i < contour.vertices.size(); ++i) {
            double dx = contour.vertices[i].x - uniqueVertices.back().x;
            double dy = contour.vertices[i].y - uniqueVertices.back().y;
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if (dist > 1e-10) {  // Only add if not too close to previous vertex
                uniqueVertices.push_back(contour.vertices[i]);
            }
        }
        
        contour.vertices = uniqueVertices;
    }
    
    std::cout << "  Total vertices: " << contour.vertices.size() << std::endl;
    
    return contour;
}

void ObstacleContour::generateDiscretization(int numSingularities) {
    (void)numSingularities; // Suppress unused parameter warning
    std::cout << "Discretizing obstacle with " << vertices.size() << " vertices..." << std::endl;
    
    singularities.clear();
    collocationPoints.clear();
    
    if (vertices.size() < 2) {
        std::cerr << "Error: Need at least 2 vertices for discretization" << std::endl;
        return;
    }
    
    // Create discrete singularities at vertices
    for (size_t i = 0; i < vertices.size(); ++i) {
        singularities.emplace_back(vertices[i], 0.0);  // Circulation will be solved for
    }
    
    // Create collocation points between singularities
    for (size_t i = 0; i < vertices.size() - 1; ++i) {
        // Midpoint between consecutive vertices (equation 7.1.24)
        double x = 0.5 * (vertices[i+1].x + vertices[i].x);
        double y = 0.5 * (vertices[i+1].y + vertices[i].y);
        Point2D collocationPos(x, y);
        
        // Calculate normal vector (outward pointing)
        double dx = vertices[i+1].x - vertices[i].x;
        double dy = vertices[i+1].y - vertices[i].y;
        double length = std::sqrt(dx*dx + dy*dy);
        
        // Normal vector: rotate tangent by 90 degrees and normalize
        Point2D normal(-dy/length, dx/length);
        
        collocationPoints.emplace_back(collocationPos, normal);
    }
    
    std::cout << "Created " << singularities.size() << " singularities and " 
              << collocationPoints.size() << " collocation points" << std::endl;
}



// Linear System Solver Implementation

LinearSystemSolver::LinearSystemSolver(int numSingularities) : M(numSingularities) {
    // System size: M-1 equations from non-penetration + 1 equation from circulation
    A = Eigen::MatrixXd::Zero(M, M);
    b = Eigen::VectorXd::Zero(M);
}

void LinearSystemSolver::buildSystem(const std::vector<DiscreteSingularity>& singularities,
                                   const std::vector<CollocationPoint>& collocationPoints,
                                   const Point2D& freeStreamVelocity,
                                   double totalCirculation) {
    
    std::cout << "Building " << M << "x" << M << " linear system..." << std::endl;
    
    // Equation 7.1.19: Non-penetration condition at collocation points
    for (size_t k = 0; k < collocationPoints.size(); ++k) {
        const auto& colloc = collocationPoints[k];
        
        for (size_t j = 0; j < singularities.size(); ++j) {
            const auto& singularity = singularities[j];
            
            // Calculate velocity induced by singularity j at collocation point k
            Point2D inducedVel = vortexVelocity(
                colloc.position.x, colloc.position.y,
                singularity.position.x, singularity.position.y,
                1.0  // Unit circulation for matrix element
            );
            
            // Matrix element: (n_k, V_j)
            double matrixElement = colloc.normal.x * inducedVel.x + colloc.normal.y * inducedVel.y;
            A(k, j) = matrixElement;
        }
        
        // Right-hand side: -n_k · V_∞
        b(k) = -(colloc.normal.x * freeStreamVelocity.x + colloc.normal.y * freeStreamVelocity.y);
    }
    
    // Equation 7.1.20: Circulation preservation condition
    for (size_t j = 0; j < singularities.size(); ++j) {
        A(collocationPoints.size(), j) = 1.0;  // Sum of circulations
    }
    b(collocationPoints.size()) = totalCirculation;
    
    // Check matrix condition
    double det = A.determinant();
    if (std::abs(det) < 1e-12) {
        std::cout << "Warning: Matrix is nearly singular (det=" << det << ")" << std::endl;
    }
}

std::vector<double> LinearSystemSolver::solve() {
    std::cout << "Solving with QR decomposition..." << std::endl;
    
    // Solve the linear system AΓ = b
    Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);
    
    // Convert to std::vector
    std::vector<double> circulations;
    circulations.reserve(M);
    for (int i = 0; i < M; ++i) {
        circulations.push_back(solution(i));
    }
    
    // Check solution quality
    Eigen::VectorXd residual = A * solution - b;
    double relativeError = residual.norm() / b.norm();
    
    if (relativeError > 1e-6) {
        std::cout << "Warning: High solution error (" << relativeError << ")" << std::endl;
    }
    
    return circulations;
}

} // namespace CFD 