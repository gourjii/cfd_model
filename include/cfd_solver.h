#pragma once

#include "discrete_singularity.h"
#include <vector>

namespace CFD {

// Field data structure for visualization
struct FieldData {
    std::vector<std::vector<double>> u;  // x-velocity component
    std::vector<std::vector<double>> v;  // y-velocity component
    std::vector<std::vector<double>> potential;  // velocity potential φ
    std::vector<std::vector<double>> streamFunction;  // stream function ψ
    std::vector<std::vector<double>> velocityMagnitude;  // |V|
    std::vector<std::vector<double>> vorticity;  // vorticity
    int nx, ny;
    double dx, dy;
    
    FieldData(int nx, int ny, double dx, double dy);
};

// Main CFD solver class
class CFDSolver {
protected:
    ObstacleContour obstacle;
    std::vector<DiscreteSingularity> singularities;
    std::vector<CollocationPoint> collocationPoints;
    Point2D freeStreamVelocity;
    double totalCirculation;
    FieldData fieldData;
    
public:
    CFDSolver(int nx, int ny, double dx, double dy);
    
    // Set up the problem
    void setupFromLineSegments(const std::vector<std::pair<Point2D, Point2D>>& lineSegments, int numSingularities);
    void setFreeStreamVelocity(double u, double v);
    void setTotalCirculation(double circulation);
    
    // Solve the problem
    virtual void solve();
    
    // Calculate field data
    virtual void calculateFieldData();
    
    // Get results
    const std::vector<DiscreteSingularity>& getSingularities() const { return singularities; }
    const ObstacleContour& getObstacleContour() const { return obstacle; }
    
    // Get flow parameters
    double getFreeStreamU() const { return freeStreamVelocity.x; }
    double getFreeStreamV() const { return freeStreamVelocity.y; }
    double getTotalCirculation() const { return totalCirculation; }
    
    // Update boundary vortex circulations (for unsteady solver)
    void updateBoundaryVortexCirculations(const std::vector<double>& circulations);
    
    // Write VTK file for visualization
    void writeVTKFile(const std::string& filename);
    void writeVTKFile(const std::string& filename, double timeValue);
    
    // Write obstacle as thick tube representation for better visibility
    void writeThickObstacleVTKFile(const std::string& filename);
    
    // Debug information
    void printSolutionInfo() const;
};

} // namespace CFD 