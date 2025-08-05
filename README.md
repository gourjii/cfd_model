# CFD Potential Flow Solver

A C++ implementation of the **Method of Discrete Singularities (МДО)** for solving 2D potential flow problems around obstacles.

## Prerequisites

- LLVM/Clang compiler
- CMake (version 3.16 or higher)

## Building the Project

### Quick Build
```bash
./build.sh
```

### Manual Build
```bash
mkdir build
cd build
cmake ..
make
```

## Running the CFD Solver

After building, you can run:

```bash
./build/bin/potential_flow_solver
```

This implements the complete CFD solution using the **Method of Discrete Singularities (МДО)** as described in the Ukrainian CFD problem. It solves the potential flow around a vertical line obstacle using:
- **Linear system solution** (equations 7.1.19 and 7.1.20)
- **Complex potential theory** with regularization
- **Multiple field outputs**: velocity, potential, stream function, vorticity

### Viewing the Results
```bash
./view_potential_flow.sh
```

This opens ParaView with the potential flow solution using the discrete singularity method.

For complete visualization including the obstacle geometry:
```bash
./view_with_obstacle.sh
```

This opens ParaView with both the flow field and obstacle geometry files.

## Project Structure

```
cfd_model/
├── CMakeLists.txt              # CMake configuration
├── src/                        # Source files
│   ├── potential_flow_solver.cpp # Main discrete singularity solver
│   ├── cfd_solver.cpp          # Core CFD solver implementation
│   ├── discrete_singularity.cpp # Discrete singularity method
│   └── complex_math.cpp        # Mathematical utilities
├── include/                    # Header files
│   ├── cfd_solver.h
│   ├── complex_math.h
│   └── discrete_singularity.h
├── build.sh                   # Build script
├── view_potential_flow.sh     # ParaView launcher (flow field only)
├── view_with_obstacle.sh      # ParaView launcher (flow + obstacle)
├── output/                    # VTK output files with timestamps
│   ├── YYYYMMDD_HHMMSS_potential_flow_discrete_singularity.vtk  # Flow field data
│   └── YYYYMMDD_HHMMSS_obstacle_geometry.vtk                    # Obstacle geometry
└── README.md                  # This file
```

## Mathematical Foundation

This implementation solves the 2D potential flow problem using:

### Governing Equations
- **Laplace equation**: Δφ = 0 (incompressible, irrotational flow)
- **Boundary conditions**: Non-penetration on obstacle surface
- **Discrete singularity method**: Point vortices along obstacle contour

### Key Components
- **Complex potential theory** with regularization
- **Linear system solver** for circulation strengths
- **Field reconstruction** across entire domain
- **Multiple visualization outputs**: velocity, potential, streamlines, vorticity
- **Obstacle visualization**: Separate VTK file for obstacle geometry

## Development

This project uses:
- C++17 standard
- CMake for build configuration
- LLVM/Clang compiler
- Eigen3 for linear algebra
- Warning flags for better code quality