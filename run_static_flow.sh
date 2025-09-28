#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PROJECT_ROOT="/Users/roman.hurzhii/WebstormProjects/cfd_model"
BUILD_DIR="$PROJECT_ROOT/build"
OUTPUT_DIR="$PROJECT_ROOT/output/static_flow"
VENV_DIR="$PROJECT_ROOT/cfd_venv"

echo -e "${BLUE}=== CFD Static Flow Simulation ===${NC}"
echo -e "${BLUE}Project: $PROJECT_ROOT${NC}"

# Step 1: Build the project
echo -e "\n${YELLOW}Step 1: Building CFD project...${NC}"
cd "$BUILD_DIR" || exit 1
make
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Build successful${NC}"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi

# Step 2: Clean old output
echo -e "\n${YELLOW}Step 2: Cleaning old static output...${NC}"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
echo -e "${GREEN}✓ Output directory cleaned${NC}"

# Step 3: Run static simulation
echo -e "\n${YELLOW}Step 3: Running static CFD simulation...${NC}"
cd "$BUILD_DIR"
./bin/potential_flow_solver

# Count generated files
VTK_COUNT=$(find "$OUTPUT_DIR" -name "*.vtk" -type f | wc -l | tr -d ' ')
echo -e "${GREEN}✓ Static simulation completed - Generated $VTK_COUNT VTK files${NC}"

# Step 4: Create visualization
echo -e "\n${YELLOW}Step 4: Creating static visualization...${NC}"

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Find the latest VTK file
LATEST_VTK=$(find "$OUTPUT_DIR" -name "*potential_flow*.vtk" -type f | sort | tail -1)

if [[ -n "$LATEST_VTK" ]]; then
    echo -e "${BLUE}Creating static visualization...${NC}"
    
    # Use the main visualization script
    cd "$PROJECT_ROOT"
    python3 - << EOF
import sys
sys.path.append('.')
from matplotlib_cfd_viewer import parse_vtk_structured_grid, extract_obstacle_from_vtk_data, create_2d_visualization

# Read the VTK file
vtk_file = "$LATEST_VTK"
output_path = "$OUTPUT_DIR/static_flow_analysis.png"

data = parse_vtk_structured_grid(vtk_file)
obstacle_data = extract_obstacle_from_vtk_data(data)

# Create visualization
fig = create_2d_visualization(data, 'CFD Static Flow Analysis', output_path, obstacle_data)
print(f"Static visualization saved to: {output_path}")
EOF
    
    echo -e "${GREEN}✓ Static visualization created: $OUTPUT_DIR/static_flow_analysis.png${NC}"
else
    echo -e "${RED}✗ No VTK files found for visualization${NC}"
fi

# Summary
echo -e "\n${GREEN}=== CFD Static Simulation Complete! ===${NC}"
echo -e "${YELLOW}Generated Files:${NC}"
echo "  • VTK files: $OUTPUT_DIR/*.vtk"
echo "  • Visualization: $OUTPUT_DIR/static_flow_analysis.png"

echo -e "\n${YELLOW}Quick Commands:${NC}"
echo -e "  • Re-run simulation: ${BLUE}./run_static_flow.sh${NC}"
echo -e "  • View unsteady flow: ${BLUE}./run_unsteady_flow.sh${NC}"

echo -e "\n${GREEN}🎉 Static flow analysis complete!${NC}"