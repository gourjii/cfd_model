#!/bin/bash

# CFD Static Flow Simulation with Matplotlib Visualization
# This script builds, runs the static simulation, and creates visualizations

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PROJECT_DIR="/Users/roman.hurzhii/WebstormProjects/cfd_model"
BUILD_DIR="$PROJECT_DIR/build"
OUTPUT_DIR="$PROJECT_DIR/output/static_flow"
VENV_DIR="$PROJECT_DIR/cfd_venv"

echo -e "${BLUE}=== CFD Static Flow Simulation ===${NC}"
echo -e "${BLUE}Project: $PROJECT_DIR${NC}"

# Step 1: Build the project
echo -e "\n${YELLOW}Step 1: Building CFD project...${NC}"
cd "$BUILD_DIR"

if make; then
    echo -e "${GREEN}✓ Build successful${NC}"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi

# Step 2: Clean old output
echo -e "\n${YELLOW}Step 2: Cleaning old static output...${NC}"
rm -f "$OUTPUT_DIR"/*.vtk
rm -f "$OUTPUT_DIR"/*.png
echo -e "${GREEN}✓ Output directory cleaned${NC}"

# Step 3: Run static simulation
echo -e "\n${YELLOW}Step 3: Running static CFD simulation...${NC}"
if ./bin/potential_flow_solver; then
    VTK_COUNT=$(find "$OUTPUT_DIR" -name "*.vtk" | wc -l)
    echo -e "${GREEN}✓ Static simulation completed - Generated $VTK_COUNT VTK files${NC}"
else
    echo -e "${RED}✗ Static simulation failed${NC}"
    exit 1
fi

# Step 4: Create matplotlib visualization
echo -e "\n${YELLOW}Step 4: Creating static visualization...${NC}"
cd "$PROJECT_DIR"

# Activate virtual environment and create static visualization
source "$VENV_DIR/bin/activate"

# Create a simple static viewer for the flow field VTK file (not obstacle geometry)
LATEST_VTK=$(find "$OUTPUT_DIR" -name "*potential_flow*.vtk" -type f | sort | tail -1)

if [[ -n "$LATEST_VTK" ]]; then
    echo -e "${BLUE}Creating static visualization...${NC}"
    
    # Create a static visualization script
    python3 - << EOF
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import sys

def parse_vtk_structured_grid(filename):
    """Parse VTK structured grid file manually"""
    print(f"Reading: {filename}")
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    dimensions = None
    points_start = None
    
    for i, line in enumerate(lines):
        if line.startswith('DIMENSIONS'):
            dims = [int(x) for x in line.split()[1:]]
            dimensions = dims
            print(f"Dimensions: {dimensions}")
        elif line.startswith('POINTS'):
            points_start = i + 1
            num_points = int(line.split()[1])
            print(f"Number of points: {num_points}")
            break
    
    if dimensions is None or points_start is None:
        print("Could not parse VTK file structure")
        return None
    
    # Read points
    points = []
    for i in range(points_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('FIELD') or line.startswith('POINT_DATA') or line.startswith('VECTORS'):
            break
        coords = [float(x) for x in line.split()]
        for j in range(0, len(coords), 3):
            if j + 2 < len(coords):
                points.append([coords[j], coords[j+1], coords[j+2]])
    
    points = np.array(points)
    print(f"Read {len(points)} points")
    
    # Read velocity data
    velocity_data = []
    velocity_start = None
    
    for i, line in enumerate(lines):
        if 'VECTORS velocity' in line:
            velocity_start = i + 1
            break
    
    if velocity_start is not None:
        for i in range(velocity_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('SCALARS') or line.startswith('LOOKUP_TABLE'):
                break
            velocities = [float(x) for x in line.split()]
            for j in range(0, len(velocities), 3):
                if j + 2 < len(velocities):
                    velocity_data.append([velocities[j], velocities[j+1], velocities[j+2]])
    
    velocity_data = np.array(velocity_data) if velocity_data else None
    
    if velocity_data is not None:
        print(f"Read {len(velocity_data)} velocity vectors")
        velocity_magnitude = np.sqrt(np.sum(velocity_data**2, axis=1))
    else:
        print("No velocity data found")
        velocity_magnitude = None
    
    return {
        'points': points,
        'dimensions': dimensions,
        'velocity': velocity_data,
        'velocity_magnitude': velocity_magnitude
    }

def create_static_visualization(vtk_file, output_path):
    """Create static visualization of CFD data"""
    data = parse_vtk_structured_grid(vtk_file)
    if data is None or data['points'] is None:
        return False
    
    points = data['points']
    dims = data['dimensions']
    velocity = data['velocity']
    velocity_mag = data['velocity_magnitude']
    
    # Extract 2D coordinates
    x = points[:, 0].reshape(dims[1], dims[0])
    y = points[:, 1].reshape(dims[1], dims[0])
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('CFD Static Flow Analysis', fontsize=16)
    
    # Plot 1: Velocity magnitude contour
    if velocity_mag is not None:
        vel_mag_2d = velocity_mag.reshape(dims[1], dims[0])
        
        # Contour plot
        contour = ax1.contourf(x, y, vel_mag_2d, levels=20, cmap='jet')
        ax1.set_title('Velocity Magnitude')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_aspect('equal')
        plt.colorbar(contour, ax=ax1, label='Velocity Magnitude')
        
        # Add contour lines
        ax1.contour(x, y, vel_mag_2d, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    
    # Plot 2: Velocity vectors
    if velocity is not None:
        u = velocity[:, 0].reshape(dims[1], dims[0])
        v = velocity[:, 1].reshape(dims[1], dims[0])
        
        # Subsample for vector plot
        step = max(1, dims[0] // 20)
        x_sub = x[::step, ::step]
        y_sub = y[::step, ::step]
        u_sub = u[::step, ::step]
        v_sub = v[::step, ::step]
        
        # Vector plot with velocity magnitude as background
        if velocity_mag is not None:
            vel_mag_2d = velocity_mag.reshape(dims[1], dims[0])
            im = ax2.imshow(vel_mag_2d, extent=[x.min(), x.max(), y.min(), y.max()], 
                           origin='lower', cmap='jet', alpha=0.7)
            plt.colorbar(im, ax=ax2, label='Velocity Magnitude')
        
        # Add velocity vectors
        ax2.quiver(x_sub, y_sub, u_sub, v_sub, scale=None, alpha=0.8, color='white', width=0.003)
        ax2.set_title('Velocity Vectors')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.show()
    
    return True

# Main execution
vtk_file = "$LATEST_VTK"
output_path = "$OUTPUT_DIR/static_flow_analysis.png"
create_static_visualization(vtk_file, output_path)
EOF

    echo -e "${GREEN}✓ Static visualization created: $OUTPUT_DIR/static_flow_analysis.png${NC}"
else
    echo -e "${RED}✗ No VTK files found for visualization${NC}"
fi

# Step 5: Summary
echo -e "\n${GREEN}=== CFD Static Simulation Complete! ===${NC}"
echo -e "${YELLOW}Generated Files:${NC}"
echo -e "  • VTK files: $OUTPUT_DIR/*.vtk"
echo -e "  • Visualization: $OUTPUT_DIR/static_flow_analysis.png"
echo -e ""
echo -e "${YELLOW}Quick Commands:${NC}"
echo -e "  • Re-run simulation: ${BLUE}./run_static_flow.sh${NC}"
echo -e "  • View unsteady flow: ${BLUE}./run_unsteady_flow.sh${NC}"
echo -e ""
echo -e "${GREEN}🎉 Static flow analysis complete!${NC}"
