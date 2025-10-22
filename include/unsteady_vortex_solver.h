#pragma once

#include "cfd_solver.h"
#include "moving_vortex.h"
#include "discrete_singularity.h"
#include <vector>
#include <memory>
#include <string>

namespace CFD {

// ========================================================================
// Extended linear system solver for time-dependent problems
// Implements equation (31): System of equations for circulation strengths Γ_j(t_{n+1})
// ========================================================================
class UnsteadyLinearSystemSolver : public LinearSystemSolver {
private:
    std::vector<double> previousCirculations;  // Γ_j(t_n) - previous time step circulations
    double timeStep;                           // Δt - current time step
    
public:
    UnsteadyLinearSystemSolver(int numSingularities);
    
    // Build system with time-dependent terms (Equation 31)
    // Solves: ∑_{j=1}^M Γ_j(t_{n+1})[V̄_j · n̄_k] = -[V̄_∞ · n̄_k] - ∑_p ∑_i γ_i^p[V̄_i^p · n̄_k]
    void buildUnsteadySystem(const std::vector<DiscreteSingularity>& boundaryVortices,
                           const std::vector<MovingVortex>& wakeVortices,
                           const std::vector<CollocationPoint>& collocationPoints,
                           const Point2D& freeStreamVelocity,
                           double totalCirculation,
                           double currentTime,
                           double timeStep);
    
    // Set previous time step circulations
    void setPreviousCirculations(const std::vector<double>& prevCirc) { 
        previousCirculations = prevCirc; 
    }
    
    // Set time step
    void setTimeStep(double dt) { timeStep = dt; }
};

// Main unsteady vortex solver class
class UnsteadyVortexSolver : public CFDSolver {
private:
    // Time integration
    double currentTime;
    double startTime;
    TimeIntegrationParams timeParams;
    
    // Vortex management
    std::unique_ptr<VortexTracker> vortexTracker;
    VortexSheddingParams sheddingParams;
    double lastSheddingTime;
    
    // Dynamic solver
    std::unique_ptr<UnsteadyLinearSystemSolver> unsteadySolver;
    std::vector<double> previousCirculations;
    
    // Output management
    int outputCounter;
    std::string baseOutputName;
    
    // Helper methods for vortex dynamics
    void updateVortexPositions();
    void shedNewVortices();
    bool shouldShedVortex() const;
    Point2D calculateVortexVelocity(const Point2D& position) const;
    Point2D calculateSheddingLocation() const;
    double calculateSheddingStrength() const;
    
    // Helper methods for time integration
    double calculateAdaptiveTimeStep() const;
    void updateBoundaryConditions();
    
    // Output methods
    void writeTimeStepOutput();
    std::string generateTimestampedFilename(const std::string& suffix) const;
    
public:
    // Constructor
    UnsteadyVortexSolver(int nx, int ny, double dx, double dy);
    
    // Destructor
    virtual ~UnsteadyVortexSolver() = default;
    
    // Setup methods
    void setTimeIntegrationParams(const TimeIntegrationParams& params);
    void setVortexSheddingParams(const VortexSheddingParams& params);
    void setBaseOutputName(const std::string& name) { baseOutputName = name; }
    
    // Main simulation methods
    void initializeUnsteadySolution();
    void runUnsteadySimulation();
    void stepForward(); // Single time step
    
    // Override solve method for unsteady case
    void solve() override;
    
    // Getters for analysis
    double getCurrentTime() const { return currentTime; }
    const VortexTracker& getVortexTracker() const { return *vortexTracker; }
    const TimeIntegrationParams& getTimeParams() const { return timeParams; }
    
    // Analysis methods
    void calculateUnsteadyPressureField();
    double calculateTotalWakeCirculation() const;
    Point2D calculateCenterOfVorticity() const;
    
    // Override field calculation to include time-dependent terms
    void calculateFieldData() override;
    
    // VTK output with wake vortices
    void writeUnsteadyVTKFile(const std::string& filename);
    void writeWakeVorticesVTKFile(const std::string& filename);
    
    // Debug and monitoring
    void printTimeStepInfo() const;
    void printWakeStatistics() const;
};

} // namespace CFD

