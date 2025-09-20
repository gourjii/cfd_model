#include "unsteady_vortex_solver.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "CFD Unsteady Vortex Solver - Dynamic Vortex Shedding" << std::endl;
    std::cout << "=====================================================" << std::endl;
    
    // ========================================================================
    // INPUT PARAMETERS
    // ========================================================================
    std::cout << "\n=== INPUT PARAMETERS ===" << std::endl;
    
    // Grid parameters
    int nx = 100, ny = 50;
    double dx = 0.1, dy = 0.1;
    std::cout << "Grid: " << nx << "x" << ny << " points, spacing " << dx << "x" << dy << " units" << std::endl;
    
    // Create unsteady vortex solver
    CFD::UnsteadyVortexSolver solver(nx, ny, dx, dy);
    
    // Define obstacle shape (same as static case for comparison)
    std::vector<std::pair<CFD::Point2D, CFD::Point2D>> lineSegments = {
        // vertical line for testing
        {CFD::Point2D(0, 0.5), CFD::Point2D(0, -0.5)}
    };
    
    int numSingularities = 30;
    solver.setupFromLineSegments(lineSegments, numSingularities);
    
    // Flow parameters
    double magnitude = 1.0;
    double angle_degrees = 0.0;  // Horizontal flow for now
    double angle_radians = angle_degrees * M_PI / 180.0;
    double freeStreamU = magnitude * cos(angle_radians);
    double freeStreamV = magnitude * sin(angle_radians);
    double totalCirculation = 0.0;
    
    solver.setFreeStreamVelocity(freeStreamU, freeStreamV);
    solver.setTotalCirculation(totalCirculation);
    
    // Time integration parameters
    CFD::TimeIntegrationParams timeParams;
    timeParams.timeStep = 0.02;
    timeParams.totalTime = 1.0;  // Shorter test run
    timeParams.outputInterval = 0.1;
    timeParams.adaptiveTimeStep = false;
    
    solver.setTimeIntegrationParams(timeParams);
    
    // Vortex shedding parameters
    CFD::VortexSheddingParams sheddingParams;
    sheddingParams.enableShedding = true;
    sheddingParams.sheddingStrength = 0.5;
    sheddingParams.sheddingInterval = 0.2;
    sheddingParams.sheddingLocation = CFD::Point2D(0.6, 0.0);  // Behind the obstacle
    
    solver.setVortexSheddingParams(sheddingParams);
    
    // Set base output name
    solver.setBaseOutputName("unsteady_vortex_flow");
    
    std::cout << "Flow: free stream (" << freeStreamU << "," << freeStreamV << ")" << std::endl;
    std::cout << "Time integration: dt=" << timeParams.timeStep << ", total=" << timeParams.totalTime << std::endl;
    std::cout << "Vortex shedding: strength=" << sheddingParams.sheddingStrength 
              << ", interval=" << sheddingParams.sheddingInterval << std::endl;
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    
    try {
        // Initialize and run unsteady simulation
        solver.solve(); // This calls initializeUnsteadySolution()
        solver.runUnsteadySimulation();
        
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "SIMULATION COMPLETED SUCCESSFULLY" << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
        // Print final statistics
        solver.printWakeStatistics();
        
        std::cout << "\nOutput files written to output/ directory" << std::endl;
        std::cout << "Use ParaView to visualize the unsteady flow and wake vortices" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during simulation: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

