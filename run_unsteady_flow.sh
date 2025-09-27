#!/bin/bash

# CFD Unsteady Simulation with Matplotlib Visualization
# This script builds, runs the simulation, and creates visualizations

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
OUTPUT_DIR="$PROJECT_DIR/output/unsteady_vortex"
VENV_DIR="$PROJECT_DIR/cfd_venv"

echo -e "${BLUE}=== CFD Unsteady Flow Simulation ===${NC}"
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
echo -e "\n${YELLOW}Step 2: Cleaning old output...${NC}"
rm -f "$OUTPUT_DIR"/*.vtk
rm -f "$OUTPUT_DIR"/*.gif
rm -rf "$OUTPUT_DIR"/frame_images
echo -e "${GREEN}✓ Output directory cleaned${NC}"

# Step 3: Run simulation
echo -e "\n${YELLOW}Step 3: Running unsteady CFD simulation...${NC}"
if ./bin/unsteady_vortex_solver; then
    VTK_COUNT=$(find "$OUTPUT_DIR" -name "*.vtk" | wc -l)
    echo -e "${GREEN}✓ Simulation completed - Generated $VTK_COUNT VTK files${NC}"
else
    echo -e "${RED}✗ Simulation failed${NC}"
    exit 1
fi

# Step 4: Activate virtual environment and create visualizations
echo -e "\n${YELLOW}Step 4: Creating visualizations...${NC}"
cd "$PROJECT_DIR"

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Create animated GIF
echo -e "${BLUE}Creating animated GIF...${NC}"
if python matplotlib_cfd_viewer.py --animate --fps 4; then
    echo -e "${GREEN}✓ Animation created: $OUTPUT_DIR/cfd_animation.gif${NC}"
else
    echo -e "${RED}✗ Animation creation failed${NC}"
fi

# Create individual frame images
echo -e "${BLUE}Creating individual frame images...${NC}"
if python matplotlib_cfd_viewer.py --save-frames; then
    FRAME_COUNT=$(find "$OUTPUT_DIR/frame_images" -name "*.png" | wc -l)
    echo -e "${GREEN}✓ Created $FRAME_COUNT frame images in $OUTPUT_DIR/frame_images/${NC}"
else
    echo -e "${RED}✗ Frame image creation failed${NC}"
fi

# Step 5: Display results
echo -e "\n${YELLOW}Step 5: Opening results...${NC}"

# Open the animated GIF
if [[ -f "$OUTPUT_DIR/cfd_animation.gif" ]]; then
    echo -e "${BLUE}Opening animated GIF...${NC}"
    open "$OUTPUT_DIR/cfd_animation.gif" &
fi

# Show a static visualization
echo -e "${BLUE}Creating interactive static visualization...${NC}"
python matplotlib_cfd_viewer.py --frame 5 &

# Step 6: Summary
echo -e "\n${GREEN}=== CFD Simulation Complete! ===${NC}"
echo -e "${YELLOW}Generated Files:${NC}"
echo -e "  • VTK files: $OUTPUT_DIR/*.vtk"
echo -e "  • Animated GIF: $OUTPUT_DIR/cfd_animation.gif"
echo -e "  • Frame images: $OUTPUT_DIR/frame_images/*.png"
echo -e ""
echo -e "${YELLOW}Quick Commands:${NC}"
echo -e "  • Re-run unsteady: ${BLUE}./run_unsteady_flow.sh${NC}"
echo -e "  • Run static flow: ${BLUE}./run_static_flow.sh${NC}"
echo -e "  • Quick animation: ${BLUE}./quick_animate.sh${NC}"
echo -e "  • View frame: ${BLUE}./quick_view.sh [N]${NC}"
echo -e ""
echo -e "${GREEN}🎉 Ready for demo! The animated GIF should be opening now.${NC}"
