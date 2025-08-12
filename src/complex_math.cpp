#include "complex_math.h"
#include <cmath>

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



} // namespace CFD 