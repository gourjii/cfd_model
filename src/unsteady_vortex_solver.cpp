#include "unsteady_vortex_solver.h"
#include "complex_math.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <limits>

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
    
    // ========================================================================
    // FORMULA REFERENCE: Equation (31) from mathematical documentation
    // System of equations for circulation strengths Γ_j(t_{n+1}):
    // ∑_{j=1}^M Γ_j(t_{n+1})[V̄(x_k,y_k,x_{0j},y_{0j}) · n̄(x_k,y_k)] = 
    // -[V̄_∞ · n̄(x_k,y_k)] - ∑_p ∑_{i=1}^{n+1} γ_i^p[V̄(x_k,y_k,x_i^p(t_{n+1}),y_i^p(t_{n+1})) · n̄(x_k,y_k)]
    // k = 1, M-1
    // ∑_{j=1}^M Γ_j(t_{n+1}) = -∑_p ∑_{i=1}^{n+1} γ_i^p
    // ========================================================================
    
    // Non-penetration condition at collocation points (Equation 31, first part)
    for (size_t k = 0; k < collocationPoints.size(); ++k) {
        const auto& colloc = collocationPoints[k];
        
        // Left-hand side: Boundary vortex contributions ∑_{j=1}^M Γ_j(t_{n+1})[V̄_j · n̄_k]
        for (size_t j = 0; j < boundaryVortices.size(); ++j) {
            const auto& singularity = boundaryVortices[j];
            
            // Calculate velocity induced by boundary singularity j at collocation point k
            // Using equations (34-36) for velocity calculation
            Point2D inducedVel = vortexVelocity(
                colloc.position.x, colloc.position.y,
                singularity.position.x, singularity.position.y,
                1.0  // Unit circulation for matrix element
            );
            
            // Matrix element: V̄(x_k,y_k,x_{0j},y_{0j}) · n̄(x_k,y_k)
            double matrixElement = colloc.normal.x * inducedVel.x + colloc.normal.y * inducedVel.y;
            A(k, j) = matrixElement;
        }
        
        // Right-hand side: -[V̄_∞ · n̄(x_k,y_k)] - wake contributions
        double rhsElement = -(colloc.normal.x * freeStreamVelocity.x + colloc.normal.y * freeStreamVelocity.y);
        
        // Add wake vortex contributions: -∑_p ∑_{i=1}^{n+1} γ_i^p[V̄(x_k,y_k,x_i^p,y_i^p) · n̄_k]
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
    
    // ========================================================================
    // FORMULA REFERENCE: Equation (31) from mathematical documentation
    // Circulation preservation condition (second part):
    // ∑_{j=1}^M Γ_j(t_{n+1}) = -∑_p ∑_{i=1}^{n+1} γ_i^p
    // ========================================================================
    
    // Circulation preservation constraint
    for (size_t j = 0; j < boundaryVortices.size(); ++j) {
        A(collocationPoints.size(), j) = 1.0;  // Sum of boundary circulations
    }
    
    // Account for wake circulation in total circulation constraint
    double wakeCirculation = 0.0;
    for (const auto& wakeVortex : wakeVortices) {
        wakeCirculation += wakeVortex.circulation;
    }
    
    // Total circulation = boundary circulation + wake circulation = constant
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

// ========================================================================
// FORMULA REFERENCE: Complete time stepping algorithm
// Implements the full unsteady vortex method time integration from t_n to t_{n+1}
// Following the sequence: Position update (Eq. 32) → System solve (Eq. 31) → Field update
// ========================================================================
void UnsteadyVortexSolver::stepForward() {
    // Calculate adaptive time step if enabled (Equation 37)
    double dt = timeParams.adaptiveTimeStep ? calculateAdaptiveTimeStep() : timeParams.timeStep;
    
    // Step 1: Update vortex positions using equation (32): r̄_{n+1} = r̄_n + V̄_n*τ_n
    updateVortexPositions();
    
    // Step 2: Check for vortex shedding (new vortex creation)
    if (shouldShedVortex()) {
        shedNewVortices();
    }
    
    // Step 3: Solve for new circulation strengths using equation (31)
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
    
    // Solve linear system for Γ_j(t_{n+1})
    std::vector<double> newCirculations = unsteadySolver->solve();
    
    // Step 4: Update boundary vortex circulations with new values
    updateBoundaryVortexCirculations(newCirculations);
    
    // Store for next time step
    previousCirculations = newCirculations;
    unsteadySolver->setPreviousCirculations(previousCirculations);
    
    // Step 5: Advance time: t_{n+1} = t_n + τ_n
    currentTime += dt;
    vortexTracker->setCurrentTime(currentTime);
    
    // Step 6: Cleanup old or distant vortices
    Point2D domainCenter(0.0, 0.0); // Could be made configurable
    double maxDistance = 10.0; // Could be made configurable
    vortexTracker->cleanupVortices(currentTime, domainCenter, maxDistance);
    
    // Step 7: Calculate new field data using equations (33-44)
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
    // Write flow field with time information
    std::string flowFilename = generateTimestampedFilename("flow");
    writeVTKFile(flowFilename, currentTime);
    
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

// ========================================================================
// FORMULA REFERENCE: Equation (32) from mathematical documentation
// Vortex position update: r̄_{n+1} = r̄_n + V̄_n*τ_n
// where V̄_n is velocity at vortex position and τ_n is time step
// ========================================================================
void UnsteadyVortexSolver::updateVortexPositions() {
    auto& wakeVortices = vortexTracker->getWakeVorticesMutable();
    
    // Calculate adaptive time step if enabled (Equation 37)
    double dt = timeParams.adaptiveTimeStep ? calculateAdaptiveTimeStep() : timeParams.timeStep;
    
    std::cout << "Updating positions of " << wakeVortices.size() << " wake vortices..." << std::endl;
    
    for (auto& vortex : wakeVortices) {
        // Calculate velocity at vortex position using equation (33)
        Point2D velocity = calculateVortexVelocity(vortex.position);
        
        // Update position: r̄_{n+1} = r̄_n + V̄_n*τ_n (Equation 32)
        vortex.position.x += velocity.x * dt;
        vortex.position.y += velocity.y * dt;
        
        // Store velocity for potential future use
        vortex.velocity = velocity;
    }
}

// ========================================================================
// FORMULA REFERENCE: Vortex shedding mechanism
// Creates new discrete vortices in the wake to model flow separation
// Implements Karman vortex street with alternating circulation signs
// ========================================================================
void UnsteadyVortexSolver::shedNewVortices() {
    if (!sheddingParams.enableShedding) {
        return;
    }
    
    std::cout << "Shedding new vortices at t=" << currentTime << std::endl;
    
    // Get obstacle geometry to find shedding locations
    const auto& collocationPoints = getObstacleContour().getCollocationPoints();
    
    if (collocationPoints.empty()) {
        return;
    }
    
    // Calculate shedding location (typically trailing edge or corners)
    Point2D sheddingPos = calculateSheddingLocation();
    
    // Calculate circulation strength based on local flow conditions
    double circulation = calculateSheddingStrength();
    
    // Alternate sign for Karman vortex street effect
    // This creates the characteristic alternating vortex pattern
    static int sheddingCounter = 0;
    double sign = (sheddingCounter % 2 == 0) ? 1.0 : -1.0;
    sheddingCounter++;
    
    // Shed vortices from both upper and lower sides for realistic wake
    // Upper vortex (positive y offset)
    Point2D upperPos = sheddingPos;
    upperPos.y += 0.1; // Offset from centerline
    MovingVortex upperVortex(upperPos, sign * circulation, currentTime);
    vortexTracker->addVortex(upperVortex);
    
    // Lower vortex (negative y offset, opposite circulation)
    Point2D lowerPos = sheddingPos;
    lowerPos.y -= 0.1; // Offset from centerline
    MovingVortex lowerVortex(lowerPos, -sign * circulation, currentTime);
    vortexTracker->addVortex(lowerVortex);
    
    // Update last shedding time
    lastSheddingTime = currentTime;
    
    std::cout << "  Added vortex pair with circulation ±" << circulation 
              << " at (" << sheddingPos.x << ", " << sheddingPos.y << "±0.1)" << std::endl;
}

bool UnsteadyVortexSolver::shouldShedVortex() const {
    // TODO: Implement shedding criteria
    return (currentTime - lastSheddingTime) >= sheddingParams.sheddingInterval;
}

// ========================================================================
// FORMULA REFERENCE: Equations (33-36) from mathematical documentation
// Velocity calculation: V̄(x,y,t_n) = ∇φ(x,y,t_n) = (u(x,y,t_n), v(x,y,t_n))
// V̄(x,y,t_n) = V̄_∞ + ∑_{j=1}^M Γ_j(t_n)V̄(x,y,x_{0j},y_{0j}) + ∑_p ∑_{i=1}^n γ_i^p V̄(x,y,x_i^p(t_n),y_i^p(t_n))
// 
// Where V̄(x,y,x_{0i},y_{0i}) = {u(x,y,x_{0i},y_{0i}) = (1/2π) * (y_{0i}-y)/R_{0i}^2
//                               {v(x,y,x_{0i},y_{0i}) = (1/2π) * (x-x_{0i})/R_{0i}^2
// 
// R_{0i} = √[(x-x_{0i})² + (y-y_{0i})²] if distance > δ, else δ (Equation 36)
// ========================================================================
Point2D UnsteadyVortexSolver::calculateVortexVelocity(const Point2D& position) const {
    Point2D velocity(0.0, 0.0);
    
    // Add free stream velocity V̄_∞ (Equation 33)
    velocity.x += getFreeStreamU();
    velocity.y += getFreeStreamV();
    
    // Add contributions from boundary vortices: ∑_{j=1}^M Γ_j(t_n)V̄(x,y,x_{0j},y_{0j})
    const auto& boundaryVortices = getObstacleContour().getSingularities();
    for (const auto& singularity : boundaryVortices) {
        Point2D inducedVel = vortexVelocity(
            position.x, position.y,
            singularity.position.x, singularity.position.y,
            singularity.circulation
        );
        velocity.x += inducedVel.x;
        velocity.y += inducedVel.y;
    }
    
    // Add contributions from wake vortices: ∑_p ∑_{i=1}^n γ_i^p V̄(x,y,x_i^p(t_n),y_i^p(t_n))
    const auto& wakeVortices = vortexTracker->getWakeVortices();
    for (const auto& wakeVortex : wakeVortices) {
        // Skip self-interaction (vortex doesn't induce velocity on itself)
        double dx = position.x - wakeVortex.position.x;
        double dy = position.y - wakeVortex.position.y;
        double distance = std::sqrt(dx*dx + dy*dy);
        
        if (distance > 1e-8) { // Avoid singularity
            Point2D inducedVel = vortexVelocity(
                position.x, position.y,
                wakeVortex.position.x, wakeVortex.position.y,
                wakeVortex.circulation
            );
            velocity.x += inducedVel.x;
            velocity.y += inducedVel.y;
        }
    }
    
    return velocity;
}

// ========================================================================
// Calculate optimal vortex shedding location
// Typically at the trailing edge or sharp corners of the obstacle
// ========================================================================
Point2D UnsteadyVortexSolver::calculateSheddingLocation() const {
    // Find trailing edge of obstacle (rightmost point for left-to-right flow)
    const auto& collocationPoints = getObstacleContour().getCollocationPoints();
    
    if (collocationPoints.empty()) {
        return sheddingParams.sheddingLocation;
    }
    
    // Find point with maximum x-coordinate (trailing edge for horizontal flow)
    double maxX = -std::numeric_limits<double>::max();
    Point2D trailingEdge;
    
    for (const auto& point : collocationPoints) {
        if (point.position.x > maxX) {
            maxX = point.position.x;
            trailingEdge = point.position;
        }
    }
    
    // Shed slightly downstream of trailing edge to avoid singularity
    trailingEdge.x += 0.05; // Small offset downstream
    
    return trailingEdge;
}

// ========================================================================
// Calculate circulation strength for newly shed vortices
// Based on local flow velocity and characteristic length scale
// Γ ≈ U * L (dimensional analysis)
// ========================================================================
double UnsteadyVortexSolver::calculateSheddingStrength() const {
    // Method 1: Use base shedding strength from parameters
    double baseStrength = sheddingParams.sheddingStrength;
    
    // Method 2: Calculate based on local velocity (more physical)
    Point2D sheddingPos = calculateSheddingLocation();
    Point2D velocity = calculateVortexVelocity(sheddingPos);
    double velocityMag = std::sqrt(velocity.x*velocity.x + velocity.y*velocity.y);
    
    // Avoid division by zero or very small velocities
    if (velocityMag < 1e-6) {
        velocityMag = getFreeStreamU(); // Use free stream as fallback
    }
    
    // Vortex strength proportional to velocity and characteristic length
    // Γ ≈ U * L (where L is characteristic length scale)
    const auto& collocationPoints = getObstacleContour().getCollocationPoints();
    double characteristicLength = 0.1; // Default value
    
    if (collocationPoints.size() > 1) {
        // Estimate characteristic length from panel sizes
        double dx = collocationPoints[1].position.x - collocationPoints[0].position.x;
        double dy = collocationPoints[1].position.y - collocationPoints[0].position.y;
        characteristicLength = std::sqrt(dx*dx + dy*dy);
    }
    
    double calculatedStrength = velocityMag * characteristicLength;
    
    // Use weighted combination of base and calculated strength
    double finalStrength = 0.3 * baseStrength + 0.7 * calculatedStrength;
    
    // Clamp to reasonable range
    finalStrength = std::max(0.01, std::min(finalStrength, 2.0));
    
    return finalStrength;
}

// ========================================================================
// FORMULA REFERENCE: Equation (37) from mathematical documentation
// Adaptive time step calculation: τ_n = min(δ_k) / max(|V̄|)
// where δ_k is characteristic length scale and |V̄| is velocity magnitude
// ========================================================================
double UnsteadyVortexSolver::calculateAdaptiveTimeStep() const {
    if (!timeParams.adaptiveTimeStep) {
        return timeParams.timeStep;
    }
    
    // Find minimum characteristic length scale (δ_k)
    double minDelta = std::numeric_limits<double>::max();
    
    // Check distances between boundary points
    const auto& collocationPoints = getObstacleContour().getCollocationPoints();
    for (size_t i = 0; i < collocationPoints.size(); ++i) {
        for (size_t j = i + 1; j < collocationPoints.size(); ++j) {
            double dx = collocationPoints[i].position.x - collocationPoints[j].position.x;
            double dy = collocationPoints[i].position.y - collocationPoints[j].position.y;
            double distance = std::sqrt(dx*dx + dy*dy);
            minDelta = std::min(minDelta, distance);
        }
    }
    
    // Find maximum velocity magnitude (max(|V̄|))
    double maxVelocity = 0.0;
    
    // Check velocities at boundary points
    for (const auto& colloc : collocationPoints) {
        Point2D velocity = calculateVortexVelocity(colloc.position);
        double velocityMag = std::sqrt(velocity.x*velocity.x + velocity.y*velocity.y);
        maxVelocity = std::max(maxVelocity, velocityMag);
    }
    
    // Check velocities of wake vortices
    const auto& wakeVortices = vortexTracker->getWakeVortices();
    for (const auto& vortex : wakeVortices) {
        double velocityMag = std::sqrt(vortex.velocity.x*vortex.velocity.x + vortex.velocity.y*vortex.velocity.y);
        maxVelocity = std::max(maxVelocity, velocityMag);
    }
    
    // Avoid division by zero
    if (maxVelocity < 1e-12) {
        maxVelocity = 1.0; // Use default if no significant velocity
    }
    
    // Calculate adaptive time step: τ_n = min(δ_k) / max(|V̄|) (Equation 37)
    double adaptiveTimeStep = minDelta / maxVelocity;
    
    // Apply safety factor and bounds
    adaptiveTimeStep *= 0.5; // Safety factor
    adaptiveTimeStep = std::max(adaptiveTimeStep, timeParams.minTimeStep);
    adaptiveTimeStep = std::min(adaptiveTimeStep, timeParams.maxTimeStep);
    
    return adaptiveTimeStep;
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

// ========================================================================
// Write wake vortices to VTK file for visualization
// Outputs vortex positions, circulation strengths, velocities, and ages
// ========================================================================
void UnsteadyVortexSolver::writeWakeVorticesVTKFile(const std::string& filename) {
    const auto& wakeVortices = vortexTracker->getWakeVortices();
    
    if (wakeVortices.empty()) {
        // Create empty file for consistency in time series
        std::ofstream file(filename);
        if (file.is_open()) {
            file << "# vtk DataFile Version 3.0\n";
            file << "Wake Vortices at time " << currentTime << " (empty)\n";
            file << "ASCII\n";
            file << "DATASET POLYDATA\n";
            file << "POINTS 0 float\n";
            file.close();
        }
        return;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Wake Vortices at time " << std::fixed << std::setprecision(6) << currentTime << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Write points (vortex positions)
    file << "POINTS " << wakeVortices.size() << " float\n";
    for (const auto& vortex : wakeVortices) {
        file << vortex.position.x << " " << vortex.position.y << " 0.0\n";
    }
    
    // Write vertices (each vortex is a point)
    file << "\nVERTICES " << wakeVortices.size() << " " << (2 * wakeVortices.size()) << "\n";
    for (size_t i = 0; i < wakeVortices.size(); ++i) {
        file << "1 " << i << "\n";
    }
    
    // Write point data (scalar and vector fields)
    file << "\nPOINT_DATA " << wakeVortices.size() << "\n";
    
    // Circulation strength (scalar)
    file << "SCALARS circulation float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vortex : wakeVortices) {
        file << vortex.circulation << "\n";
    }
    
    // Velocity vectors
    file << "\nVECTORS velocity float\n";
    for (const auto& vortex : wakeVortices) {
        file << vortex.velocity.x << " " << vortex.velocity.y << " 0.0\n";
    }
    
    // Vortex age (time since birth)
    file << "\nSCALARS age float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vortex : wakeVortices) {
        file << (currentTime - vortex.birthTime) << "\n";
    }
    
    // Absolute circulation (for coloring)
    file << "\nSCALARS abs_circulation float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& vortex : wakeVortices) {
        file << std::abs(vortex.circulation) << "\n";
    }
    
    file.close();
    std::cout << "  Wrote " << wakeVortices.size() << " wake vortices to " << filename << std::endl;
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
