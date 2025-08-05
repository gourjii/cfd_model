#include "complex_math.h"
#include <cmath>
#include <iostream>
#include <iomanip>

namespace CFD {

double angularFunction(double x, double y, double x0, double y0) {
    // Equation 7.1.14: θ_j(x,y) = (1/(2π)) arg(z - ω_0j)
    double dx = x - x0;
    double dy = y - y0;
    
    // Handle the case when point coincides with singularity
    if (std::abs(dx) < 1e-12 && std::abs(dy) < 1e-12) {
        return 0.0;
    }
    
    // Calculate angle based on quadrant
    double angle = std::atan2(dy, dx);
    
    // Normalize to [0, 1] range (instead of [-π, π])
    if (angle < 0) {
        angle += 2.0 * M_PI;
    }
    
    return angle / (2.0 * M_PI);
}

double regularizedDistance(double x, double y, double x0, double y0, double epsilon) {
    // Equation 7.1.18': R_j = max{ε_j, √((x - x_0j)² + (y - y_0j)²)}
    double distance = std::sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
    return std::max(epsilon, distance);
}

Point2D vortexVelocity(double x, double y, double x0, double y0, double circulation, double epsilon) {
    // Equation 7.1.17': V⃗_j(x,y) = ( (y_0j - y)/(2πR_j²), (x - x_0j)/(2πR_j²) )
    double R = regularizedDistance(x, y, x0, y0, epsilon);
    double factor = circulation / (2.0 * M_PI * R * R);
    
    double vx = factor * (y0 - y);
    double vy = factor * (x - x0);
    
    return Point2D(vx, vy);
}

double vortexStreamFunction(double x, double y, double x0, double y0, double circulation, double epsilon) {
    // Stream function contribution from a point vortex
    // ψ_j = -(Γ_j/(2π)) ln(R_j)
    double R = regularizedDistance(x, y, x0, y0, epsilon);
    return -(circulation / (2.0 * M_PI)) * std::log(R);
}

// Additional utility functions for detailed analysis
void analyzeVortexField(double x, double y, double x0, double y0, double circulation, double epsilon) {
    std::cout << "\n=== VORTEX FIELD ANALYSIS ===" << std::endl;
    std::cout << "Analyzing field point (" << std::fixed << std::setprecision(6) << x << ", " << y << ")" << std::endl;
    std::cout << "Induced by singularity at (" << std::fixed << std::setprecision(6) << x0 << ", " << y0 << ")" << std::endl;
    std::cout << "Circulation: " << std::fixed << std::setprecision(6) << circulation << std::endl;
    
    // Calculate distance and angle
    double dx = x - x0;
    double dy = y - y0;
    double distance = std::sqrt(dx*dx + dy*dy);
    double angle = std::atan2(dy, dx) * 180.0 / M_PI;
    
    std::cout << "\nGeometric properties:" << std::endl;
    std::cout << "  Distance: " << std::fixed << std::setprecision(6) << distance << " units" << std::endl;
    std::cout << "  Angle: " << std::fixed << std::setprecision(6) << angle << " degrees" << std::endl;
    
    // Calculate velocity field
    Point2D velocity = vortexVelocity(x, y, x0, y0, circulation, epsilon);
    double velocityMagnitude = std::sqrt(velocity.x*velocity.x + velocity.y*velocity.y);
    double velocityAngle = std::atan2(velocity.y, velocity.x) * 180.0 / M_PI;
    
    std::cout << "\nVelocity field:" << std::endl;
    std::cout << "  Velocity: (" << std::fixed << std::setprecision(6) << velocity.x << ", " << velocity.y << ")" << std::endl;
    std::cout << "  Magnitude: " << std::fixed << std::setprecision(6) << velocityMagnitude << " units" << std::endl;
    std::cout << "  Direction: " << std::fixed << std::setprecision(6) << velocityAngle << " degrees" << std::endl;
    
    // Calculate potential and stream function
    double potential = circulation * angularFunction(x, y, x0, y0);
    double streamFunction = vortexStreamFunction(x, y, x0, y0, circulation, epsilon);
    
    std::cout << "\nPotential and stream function:" << std::endl;
    std::cout << "  Potential: " << std::fixed << std::setprecision(6) << potential << std::endl;
    std::cout << "  Stream function: " << std::fixed << std::setprecision(6) << streamFunction << std::endl;
    
    // Physical interpretation
    std::cout << "\nPhysical interpretation:" << std::endl;
    if (circulation > 0) {
        std::cout << "  ✓ Positive circulation: counterclockwise flow" << std::endl;
    } else if (circulation < 0) {
        std::cout << "  ✓ Negative circulation: clockwise flow" << std::endl;
    } else {
        std::cout << "  ✓ Zero circulation: no flow" << std::endl;
    }
    
    if (distance < epsilon) {
        std::cout << "  ⚠ Point very close to singularity - regularization applied" << std::endl;
    } else {
        std::cout << "  ✓ Point sufficiently far from singularity" << std::endl;
    }
}

} // namespace CFD 