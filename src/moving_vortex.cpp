#include "moving_vortex.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace CFD {

// VortexTracker implementation
VortexTracker::VortexTracker(double minDist, double maxVortexAge) 
    : minDistance(minDist), maxAge(maxVortexAge), currentTime(0.0) {
    wakeVortices.reserve(1000); // Reserve space for efficiency
}

void VortexTracker::addVortex(const MovingVortex& vortex) {
    // Check if too close to existing vortices
    for (const auto& existingVortex : wakeVortices) {
        double dx = vortex.position.x - existingVortex.position.x;
        double dy = vortex.position.y - existingVortex.position.y;
        double distance = std::sqrt(dx*dx + dy*dy);
        
        if (distance < minDistance) {
            // Too close to existing vortex, skip addition
            return;
        }
    }
    
    // Add the vortex
    wakeVortices.push_back(vortex);
    
    std::cout << "Added wake vortex at (" << vortex.position.x << ", " << vortex.position.y 
              << ") with circulation " << vortex.circulation << std::endl;
}

void VortexTracker::updatePositions(double timeStep) {
    // This method will be called from the main solver
    // Position updates are handled in the main simulation loop
    // to account for interaction between all vortices
    (void)timeStep; // Suppress unused parameter warning
}

void VortexTracker::cleanupVortices(double currentTime, const Point2D& domainCenter, double maxDistance) {
    auto it = wakeVortices.begin();
    int removedCount = 0;
    
    while (it != wakeVortices.end()) {
        bool shouldRemove = false;
        
        // Check age
        double age = currentTime - it->birthTime;
        if (age > maxAge) {
            shouldRemove = true;
        }
        
        // Check distance from domain center
        double dx = it->position.x - domainCenter.x;
        double dy = it->position.y - domainCenter.y;
        double distance = std::sqrt(dx*dx + dy*dy);
        if (distance > maxDistance) {
            shouldRemove = true;
        }
        
        // Check for very weak vortices
        if (std::abs(it->circulation) < 1e-8) {
            shouldRemove = true;
        }
        
        if (shouldRemove) {
            it = wakeVortices.erase(it);
            removedCount++;
        } else {
            ++it;
        }
    }
    
    if (removedCount > 0) {
        std::cout << "Removed " << removedCount << " wake vortices (age/distance/weak)" << std::endl;
    }
}

double VortexTracker::getTotalWakeCirculation() const {
    double totalCirculation = 0.0;
    for (const auto& vortex : wakeVortices) {
        totalCirculation += vortex.circulation;
    }
    return totalCirculation;
}

} // namespace CFD

