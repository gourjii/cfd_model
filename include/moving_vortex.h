#pragma once

#include "complex_math.h"
#include "discrete_singularity.h"
#include <vector>
#include <memory>

namespace CFD {

// Moving vortex structure for dynamic simulations
struct MovingVortex {
    Point2D position;      // Current position (x_i^p(t), y_i^p(t))
    Point2D velocity;      // Current velocity (u, v)
    double circulation;    // Circulation strength γ_i^p
    double birthTime;      // Time when vortex was created
    int sourcePanel;       // Index of panel that shed this vortex (-1 if not from panel)
    
    MovingVortex(const Point2D& pos, double circ, double time, int panel = -1)
        : position(pos), velocity(0.0, 0.0), circulation(circ), 
          birthTime(time), sourcePanel(panel) {}
};

// Vortex tracker for managing wake vortices
class VortexTracker {
private:
    std::vector<MovingVortex> wakeVortices;
    double minDistance;     // Minimum distance between vortices
    double maxAge;          // Maximum age before vortex removal
    double currentTime;
    
public:
    VortexTracker(double minDist = 0.01, double maxVortexAge = 100.0);
    
    // Add new vortex to wake
    void addVortex(const MovingVortex& vortex);
    
    // Update positions of all wake vortices
    void updatePositions(double timeStep);
    
    // Remove old or distant vortices
    void cleanupVortices(double currentTime, const Point2D& domainCenter, double maxDistance);
    
    // Get all wake vortices
    const std::vector<MovingVortex>& getWakeVortices() const { return wakeVortices; }
    
    // Get mutable reference for position updates
    std::vector<MovingVortex>& getWakeVorticesMutable() { return wakeVortices; }
    
    // Clear all vortices
    void clear() { wakeVortices.clear(); }
    
    // Get total number of wake vortices
    size_t getVortexCount() const { return wakeVortices.size(); }
    
    // Calculate total circulation in wake
    double getTotalWakeCirculation() const;
    
    // Set current time
    void setCurrentTime(double time) { currentTime = time; }
};

// Time integration parameters
struct TimeIntegrationParams {
    double timeStep;           // Δt
    double totalTime;          // Total simulation time
    double outputInterval;     // Time between output files
    double vortexShedPeriod;   // Period for vortex shedding
    bool adaptiveTimeStep;     // Use adaptive time stepping
    double maxTimeStep;        // Maximum allowed time step
    double minTimeStep;        // Minimum allowed time step
    
    TimeIntegrationParams() 
        : timeStep(0.01), totalTime(10.0), outputInterval(0.1),
          vortexShedPeriod(0.1), adaptiveTimeStep(false),
          maxTimeStep(0.05), minTimeStep(0.001) {}
};

// Vortex shedding parameters
struct VortexSheddingParams {
    bool enableShedding;       // Enable/disable vortex shedding
    double sheddingStrength;   // Base circulation strength for shed vortices
    double sheddingThreshold;  // Velocity threshold for shedding
    Point2D sheddingLocation;  // Location where vortices are shed (relative to obstacle)
    double sheddingInterval;   // Time interval between shedding events
    
    VortexSheddingParams()
        : enableShedding(true), sheddingStrength(0.1), 
          sheddingThreshold(0.5), sheddingLocation(0.5, 0.0),
          sheddingInterval(0.1) {}
};

} // namespace CFD

