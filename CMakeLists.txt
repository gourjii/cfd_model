cmake_minimum_required(VERSION 3.16)
project(cfd_model VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags for better debugging and warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

# Find Eigen3 for linear algebra
find_package(Eigen3 REQUIRED)

# Add executables
add_executable(potential_flow_solver src/potential_flow_solver.cpp src/cfd_solver.cpp src/discrete_singularity.cpp src/complex_math.cpp)

# Set output directory
set_target_properties(potential_flow_solver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Include directories
target_include_directories(potential_flow_solver PRIVATE include)

# Link Eigen3
target_link_libraries(potential_flow_solver Eigen3::Eigen)

# Print configuration info
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}") 