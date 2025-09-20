#include "unsteady_vortex_solver.h"
#include "complex_math.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <filesystem>

namespace CFD {

// UnsteadyLinearSystemSolver implementation
UnsteadyLinearSystemSolver::UnsteadyLinearSystemSolver(int numSingularities) 
    : LinearSystemSolver(numSingularities), timeStep(0.01) {
    previousCirculations.resize(numSingularities, 0.0);
}

void UnsteadyLinearSystemSolver::buildUnsteadySystem(
    const std::vector<DiscreteSingularity>& boundaryVortices,
    const std::vector<MovingVortex>& wakeVortices,
    const std::vector<CollocationPoint>& collocationPoints,
    const Point2D& freeStreamVelocity,
    double totalCirculation,
    double currentTime,
    double timeStep) {
    
    (void)currentTime; // Suppress unused parameter warning for now
    this->timeStep = timeStep;
    
    std::cout << "Building unsteady system with " << boundaryVortices.size() 
              << " boundary vortices and " << wakeVortices.size() << " wake vortices..." << std::endl;
    
    // Get matrix and vector references (need non-const access)
    auto& A = const_cast<Eigen::MatrixXd&>(getMatrix());
    auto& b = const_cast<Eigen::VectorXd&>(getVector());
    
    // Clear system
    A.setZero();
    b.setZero();
    
    // Build system following equations from slides 15-21
    // Non-penetration condition at collocation points (modified for wake interaction)
    for (size_t k = 0; k < collocationPoints.size(); ++k) {
        const auto& colloc = collocationPoints[k];
        
        // Boundary vortex contributions
        for (size_t j = 0; j < boundaryVortices.size(); ++j) {
            const auto& singularity = boundaryVortices[j];
            
            // Calculate velocity induced by boundary singularity j at collocation point k
            Point2D inducedVel = vortexVelocity(
                colloc.position.x, colloc.position.y,
                singularity.position.x, singularity.position.y,
                1.0  // Unit circulation for matrix element
            );
            
            // Matrix element: (n_k, V_j)
            double matrixElement = colloc.normal.x * inducedVel.x + colloc.normal.y * inducedVel.y;
            A(k, j) = matrixElement;
        }
        
        // Right-hand side: -n_k · V_∞ - n_k · V_wake
        double rhsElement = -(colloc.normal.x * freeStreamVelocity.x + colloc.normal.y * freeStreamVelocity.y);
        
        // Add wake vortex contributions to RHS
        for (const auto& wakeVortex : wakeVortices) {
            Point2D wakeInducedVel = vortexVelocity(
                colloc.position.x, colloc.position.y,
                wakeVortex.position.x, wakeVortex.position.y,
                wakeVortex.circulation
            );
            
            rhsElement -= (colloc.normal.x * wakeInducedVel.x + colloc.normal.y * wakeInducedVel.y);
        }
        
        b(k) = rhsElement;
    }
    
    // Circulation preservation condition (equation from slide 16)
    for (size_t j = 0; j < boundaryVortices.size(); ++j) {
        A(collocationPoints.size(), j) = 1.0;  // Sum of circulations
    }
    
    // Account for wake circulation in total circulation constraint
    double wakeCirculation = 0.0;
    for (const auto& wakeVortex : wakeVortices) {
        wakeCirculation += wakeVortex.circulation;
    }
    
    b(collocationPoints.size()) = totalCirculation - wakeCirculation;
    
    // Check matrix condition
    double det = A.determinant();
    if (std::abs(det) < 1e-12) {
        std::cout << "Warning: Unsteady matrix is nearly singular (det=" << det << ")" << std::endl;
    }
}

// UnsteadyVortexSolver implementation
UnsteadyVortexSolver::UnsteadyVortexSolver(int nx, int ny, double dx, double dy)
    : CFDSolver(nx, ny, dx, dy), currentTime(0.0), startTime(0.0),
      lastSheddingTime(0.0), outputCounter(0), baseOutputName("unsteady_vortex") {
    
    // Initialize vortex tracker
    vortexTracker = std::make_unique<VortexTracker>();
    
    // Set default parameters
    timeParams = TimeIntegrationParams();
    sheddingParams = VortexSheddingParams();
    
    std::cout << "UnsteadyVortexSolver initialized with " << nx << "x" << ny << " grid" << std::endl;
}

void UnsteadyVortexSolver::setTimeIntegrationParams(const TimeIntegrationParams& params) {
    timeParams = params;
    std::cout << "Time integration parameters set: dt=" << params.timeStep 
              << ", total_time=" << params.totalTime << std::endl;
}

void UnsteadyVortexSolver::setVortexSheddingParams(const VortexSheddingParams& params) {
    sheddingParams = params;
    std::cout << "Vortex shedding parameters set: strength=" << params.sheddingStrength
              << ", enabled=" << (params.enableShedding ? "true" : "false") << std::endl;
}

void UnsteadyVortexSolver::initializeUnsteadySolution() {
    std::cout << "\n=== Initializing Unsteady Solution ===" << std::endl;
    
    // Reset time
    currentTime = startTime;
    lastSheddingTime = startTime;
    outputCounter = 0;
    
    // Clear wake vortices
    vortexTracker->clear();
    vortexTracker->setCurrentTime(currentTime);
    
    // Solve initial steady-state problem
    std::cout << "Solving initial steady-state problem..." << std::endl;
    CFDSolver::solve(); // Call parent class solve method
    
    // Initialize unsteady solver with proper size
    const auto& singularities = getObstacleContour().getSingularities();
    unsteadySolver = std::make_unique<UnsteadyLinearSystemSolver>(singularities.size());
    
    // Store initial circulations
    previousCirculations.clear();
    for (const auto& sing : singularities) {
        previousCirculations.push_back(sing.circulation);
    }
    
    // Initial field calculation
    calculateFieldData();
    
    // Write initial state
    writeTimeStepOutput();
    
    std::cout << "Initial solution completed. Starting time integration..." << std::endl;
}

void UnsteadyVortexSolver::solve() {
    // For unsteady case, just initialize - user should call runUnsteadySimulation()
    initializeUnsteadySolution();
}

void UnsteadyVortexSolver::stepForward() {
    // Calculate adaptive time step if enabled
    double dt = timeParams.adaptiveTimeStep ? calculateAdaptiveTimeStep() : timeParams.timeStep;
    
    // Update vortex positions first (equation 32 from slides)
    updateVortexPositions();
    
    // Check for vortex shedding
    if (shouldShedVortex()) {
        shedNewVortices();
    }
    
    // Solve for new circulation strengths
    const auto& boundaryVortices = getObstacleContour().getSingularities();
    const auto& collocationPoints = getObstacleContour().getCollocationPoints();
    
    unsteadySolver->buildUnsteadySystem(
        boundaryVortices,
        vortexTracker->getWakeVortices(),
        collocationPoints,
        Point2D(getFreeStreamU(), getFreeStreamV()),
        getTotalCirculation(),
        currentTime,
        dt
    );
    
    // Solve system
    std::vector<double> newCirculations = unsteadySolver->solve();
    
    // Update boundary vortex circulations
    updateBoundaryVortexCirculations(newCirculations);
    
    // Store for next time step
    previousCirculations = newCirculations;
    unsteadySolver->setPreviousCirculations(previousCirculations);
    
    // Update time
    currentTime += dt;
    vortexTracker->setCurrentTime(currentTime);
    
    // Cleanup old vortices
    Point2D domainCenter(0.0, 0.0); // Could be made configurable
    double maxDistance = 10.0; // Could be made configurable
    vortexTracker->cleanupVortices(currentTime, domainCenter, maxDistance);
    
    // Calculate new field data
    calculateFieldData();
}

std::string UnsteadyVortexSolver::generateTimestampedFilename(const std::string& suffix) const {
    // Create directory if it doesn't exist (relative to build directory)
    std::filesystem::create_directories("../output/unsteady_vortex");
    
    std::ostringstream oss;
    // Use zero-padded integer format for better ParaView time series recognition
    int timeStep = static_cast<int>(std::round(currentTime / timeParams.timeStep));
    oss << "../output/unsteady_vortex/" << baseOutputName << "_" << suffix 
        << "_" << std::setfill('0') << std::setw(6) << timeStep << ".vtk";
    return oss.str();
}

void UnsteadyVortexSolver::writeTimeStepOutput() {
    // Write flow field
    std::string flowFilename = generateTimestampedFilename("flow");
    writeVTKFile(flowFilename);
    
    // Write wake vortices
    std::string wakeFilename = generateTimestampedFilename("wake");
    writeWakeVorticesVTKFile(wakeFilename);
    
    outputCounter++;
    
    std::cout << "Output written at t=" << currentTime 
              << " (files: " << flowFilename << ", " << wakeFilename << ")" << std::endl;
}

void UnsteadyVortexSolver::printTimeStepInfo() const {
    std::cout << "Time: " << std::fixed << std::setprecision(4) << currentTime 
              << ", Wake vortices: " << vortexTracker->getVortexCount()
              << ", Wake circulation: " << vortexTracker->getTotalWakeCirculation() << std::endl;
}

// Placeholder implementations for helper methods (to be implemented in next phase)
void UnsteadyVortexSolver::updateVortexPositions() {
    // TODO: Implement vortex position updates using equation 32
    // For now, just a placeholder
}

void UnsteadyVortexSolver::shedNewVortices() {
    // TODO: Implement vortex shedding logic
    // For now, just a placeholder
}

bool UnsteadyVortexSolver::shouldShedVortex() const {
    // TODO: Implement shedding criteria
    return (currentTime - lastSheddingTime) >= sheddingParams.sheddingInterval;
}

Point2D UnsteadyVortexSolver::calculateVortexVelocity(const Point2D& position) const {
    // TODO: Implement velocity calculation at arbitrary point
    (void)position;
    return Point2D(0.0, 0.0);
}

Point2D UnsteadyVortexSolver::calculateSheddingLocation() const {
    // TODO: Calculate where to shed new vortices
    return sheddingParams.sheddingLocation;
}

double UnsteadyVortexSolver::calculateSheddingStrength() const {
    // TODO: Calculate circulation strength for new vortices
    return sheddingParams.sheddingStrength;
}

double UnsteadyVortexSolver::calculateAdaptiveTimeStep() const {
    // TODO: Implement adaptive time stepping
    return timeParams.timeStep;
}

void UnsteadyVortexSolver::updateBoundaryConditions() {
    // TODO: Update boundary conditions for time-dependent case
}

void UnsteadyVortexSolver::calculateUnsteadyPressureField() {
    // TODO: Implement pressure calculation with time derivatives
}

void UnsteadyVortexSolver::calculateFieldData() {
    // For now, call parent method
    // TODO: Override with unsteady pressure calculation
    CFDSolver::calculateFieldData();
}

void UnsteadyVortexSolver::writeUnsteadyVTKFile(const std::string& filename) {
    // For now, call parent method
    // TODO: Add time-dependent fields
    writeVTKFile(filename);
}

void UnsteadyVortexSolver::writeWakeVorticesVTKFile(const std::string& filename) {
    // TODO: Implement wake vortex visualization
    (void)filename;
    std::cout << "Wake vortex VTK output not yet implemented" << std::endl;
}

double UnsteadyVortexSolver::calculateTotalWakeCirculation() const {
    return vortexTracker->getTotalWakeCirculation();
}

Point2D UnsteadyVortexSolver::calculateCenterOfVorticity() const {
    // TODO: Implement center of vorticity calculation
    return Point2D(0.0, 0.0);
}

void UnsteadyVortexSolver::printWakeStatistics() const {
    std::cout << "Wake Statistics:" << std::endl;
    std::cout << "  Number of vortices: " << vortexTracker->getVortexCount() << std::endl;
    std::cout << "  Total circulation: " << vortexTracker->getTotalWakeCirculation() << std::endl;
}

void UnsteadyVortexSolver::runUnsteadySimulation() {
    std::cout << "\n=== Starting Unsteady Simulation ===" << std::endl;
    std::cout << "Total time: " << timeParams.totalTime << std::endl;
    std::cout << "Time step: " << timeParams.timeStep << std::endl;
    std::cout << "Output interval: " << timeParams.outputInterval << std::endl;
    
    double nextOutputTime = currentTime + timeParams.outputInterval;
    
    while (currentTime < timeParams.totalTime) {
        // Take time step
        stepForward();
        
        // Print progress
        if (static_cast<int>(currentTime / timeParams.outputInterval) % 10 == 0) {
            printTimeStepInfo();
        }
        
        // Write output if needed
        if (currentTime >= nextOutputTime) {
            writeTimeStepOutput();
            nextOutputTime += timeParams.outputInterval;
        }
    }
    
    // Final output
    writeTimeStepOutput();
    printWakeStatistics();
    
    std::cout << "Unsteady simulation completed!" << std::endl;
}

} // namespace CFD
