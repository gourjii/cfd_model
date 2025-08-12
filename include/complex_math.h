#pragma once

#include <complex>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace CFD {

// Complex number type
using Complex = std::complex<double>;

// 2D Point structure
struct Point2D {
    double x, y;
    Point2D(double x = 0, double y = 0) : x(x), y(y) {}
    
    // Convert to complex number
    Complex toComplex() const { return Complex(x, y); }
    
    // Distance to another point
    double distanceTo(const Point2D& other) const {
        return std::sqrt((x - other.x)*(x - other.x) + (y - other.y)*(y - other.y));
    }
};

// Angular function θ_j(x,y) from equation 7.1.14
double angularFunction(double x, double y, double x0, double y0);

// Regularized distance R_j from equation 7.1.18'
double regularizedDistance(double x, double y, double x0, double y0, double epsilon = 1e-6);

// Velocity induced by a point vortex (equation 7.1.17')
Point2D vortexVelocity(double x, double y, double x0, double y0, double circulation, double epsilon = 1e-6);

// Stream function contribution from a point vortex
double vortexStreamFunction(double x, double y, double x0, double y0, double circulation, double epsilon = 1e-6);



} // namespace CFD 