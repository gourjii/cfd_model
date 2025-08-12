#pragma once

#include "complex_math.h"
#include <vector>
#include <Eigen/Dense>

namespace CFD {

// Discrete singularity point (point vortex)
struct DiscreteSingularity {
    Point2D position;  // ω_0j = x_0j + iy_0j
    double circulation; // Γ_j
    
    DiscreteSingularity(const Point2D& pos, double circ = 0.0) 
        : position(pos), circulation(circ) {}
};

// Collocation point for boundary conditions
struct CollocationPoint {
    Point2D position;  // ω_k = x_k + iy_k
    Point2D normal;    // n_k = (n_xk, n_yk)
    
    CollocationPoint(const Point2D& pos, const Point2D& norm) 
        : position(pos), normal(norm) {}
};

// Obstacle contour defined by discrete points
class ObstacleContour {
private:
    std::vector<Point2D> vertices;
    std::vector<DiscreteSingularity> singularities;
    std::vector<CollocationPoint> collocationPoints;
    
public:
    // Create obstacle from array of line segments
    static ObstacleContour createFromLineSegments(const std::vector<std::pair<Point2D, Point2D>>& lineSegments, int numPoints);
    
    // Generate discrete singularities and collocation points
    void generateDiscretization(int numSingularities);
    
    // Get singularities
    const std::vector<DiscreteSingularity>& getSingularities() const { return singularities; }
    
    // Get collocation points
    const std::vector<CollocationPoint>& getCollocationPoints() const { return collocationPoints; }
    
    // Get obstacle vertices
    const std::vector<Point2D>& getVertices() const { return vertices; }
    

};

// Linear system solver for equations 7.1.19 and 7.1.20
class LinearSystemSolver {
private:
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    int M;  // number of singularities
    
public:
    LinearSystemSolver(int numSingularities);
    
    // Build the linear system from equations 7.1.19 and 7.1.20
    void buildSystem(const std::vector<DiscreteSingularity>& singularities,
                    const std::vector<CollocationPoint>& collocationPoints,
                    const Point2D& freeStreamVelocity,
                    double totalCirculation);
    
    // Solve the system and return circulations
    std::vector<double> solve();
    
    // Get the system matrix and vector (for debugging)
    const Eigen::MatrixXd& getMatrix() const { return A; }
    const Eigen::VectorXd& getVector() const { return b; }
};

} // namespace CFD 