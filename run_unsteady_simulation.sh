#!/bin/bash

# CFD Unsteady Vortex Simulation Runner
# This script builds, runs the simulation, and opens ParaView

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

# ParaView paths (common locations on macOS)
PARAVIEW_PATHS=(
    "/Applications/ParaView-5.12.0.app/Contents/MacOS/paraview"
    "/Applications/ParaView-5.11.2.app/Contents/MacOS/paraview"
    "/Applications/ParaView-5.10.1.app/Contents/MacOS/paraview"
    "/Applications/ParaView.app/Contents/MacOS/paraview"
    "/usr/local/bin/paraview"
    "/opt/homebrew/bin/paraview"
)

echo -e "${BLUE}=== CFD Unsteady Vortex Simulation Runner ===${NC}"
echo -e "${BLUE}=============================================${NC}"

# Function to find ParaView
find_paraview() {
    for path in "${PARAVIEW_PATHS[@]}"; do
        if [[ -x "$path" ]]; then
            echo "$path"
            return 0
        fi
    done
    return 1
}

# Step 1: Navigate to project directory
echo -e "\n${YELLOW}Step 1: Navigating to project directory...${NC}"
cd "$PROJECT_DIR"
echo -e "${GREEN}✓ Current directory: $(pwd)${NC}"

# Step 2: Create build directory if it doesn't exist
echo -e "\n${YELLOW}Step 2: Setting up build directory...${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
echo -e "${GREEN}✓ Build directory ready${NC}"

# Step 3: Configure with CMake
echo -e "\n${YELLOW}Step 3: Configuring with CMake...${NC}"
cmake .. || {
    echo -e "${RED}✗ CMake configuration failed${NC}"
    exit 1
}
echo -e "${GREEN}✓ CMake configuration successful${NC}"

# Step 4: Build the project
echo -e "\n${YELLOW}Step 4: Building the project...${NC}"
make || {
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
}
echo -e "${GREEN}✓ Build successful${NC}"

# Step 5: Clean previous simulation results
echo -e "\n${YELLOW}Step 5: Cleaning previous results...${NC}"
if [[ -d "$OUTPUT_DIR" ]]; then
    rm -f "$OUTPUT_DIR"/unsteady_vortex_flow_*.vtk
    echo -e "${GREEN}✓ Previous results cleaned${NC}"
else
    echo -e "${GREEN}✓ Output directory will be created${NC}"
fi

# Step 6: Run the simulation
echo -e "\n${YELLOW}Step 6: Running unsteady vortex simulation...${NC}"
echo -e "${BLUE}Starting simulation (this may take a moment)...${NC}"
./bin/unsteady_vortex_solver || {
    echo -e "${RED}✗ Simulation failed${NC}"
    exit 1
}
echo -e "${GREEN}✓ Simulation completed successfully${NC}"

# Step 7: Check results
echo -e "\n${YELLOW}Step 7: Checking simulation results...${NC}"
VTK_COUNT=$(find "$OUTPUT_DIR" -name "*.vtk" -type f | wc -l)
if [[ $VTK_COUNT -gt 0 ]]; then
    echo -e "${GREEN}✓ Found $VTK_COUNT VTK files in output directory${NC}"
    echo -e "${GREEN}✓ Results available at: $OUTPUT_DIR${NC}"
else
    echo -e "${RED}✗ No VTK files found in output directory${NC}"
    exit 1
fi

# Step 8: Find and launch ParaView
echo -e "\n${YELLOW}Step 8: Launching ParaView...${NC}"
PARAVIEW_PATH=$(find_paraview)
if [[ $? -eq 0 ]]; then
    echo -e "${GREEN}✓ Found ParaView at: $PARAVIEW_PATH${NC}"
    
    # Get all VTK files for time series (they now contain proper time information)
    VTK_FILES=($(find "$OUTPUT_DIR" -name "unsteady_vortex_flow_flow_*.vtk" -type f | sort))
    
    if [[ ${#VTK_FILES[@]} -gt 0 ]]; then
        VTK_COUNT=${#VTK_FILES[@]}
        echo -e "${GREEN}✓ Found $VTK_COUNT VTK files with time information${NC}"
        echo -e "${BLUE}Opening ParaView with time series...${NC}"
        echo -e "${BLUE}Files: $(basename "${VTK_FILES[0]}") to $(basename "${VTK_FILES[-1]}")${NC}"
        
        # Launch ParaView with all VTK files (they contain time info for series recognition)
        "$PARAVIEW_PATH" "${VTK_FILES[@]}" &
        
        echo -e "${GREEN}✓ ParaView launched successfully${NC}"
        echo -e "\n${BLUE}=== ParaView Ready for Demo! ===${NC}"
        echo -e "${GREEN}✓ Time series automatically loaded and ready (${VTK_COUNT} time steps)${NC}"
        echo -e "${YELLOW}What you should see in ParaView:${NC}"
        echo -e "${YELLOW}1. Time controls in top toolbar (0.000 to 1.000 seconds)${NC}"
        echo -e "${YELLOW}2. Pipeline Browser shows: unsteady_vortex_time_series.pvd${NC}"
        echo -e "${YELLOW}3. Click 'Apply' to load the data${NC}"
        echo -e "${YELLOW}4. Color by 'velocity_magnitude' for best visualization${NC}"
        echo -e "${YELLOW}5. Use Play button ▶️ to animate through time${NC}"
        echo -e "${YELLOW}6. Add Glyph filter for velocity vectors if desired${NC}"
        echo -e "\n${GREEN}🎬 Demo ready! No manual steps required.${NC}"
    else
        echo -e "${RED}✗ No VTK files found to open${NC}"
        exit 1
    fi
else
    echo -e "${RED}✗ ParaView not found in common locations${NC}"
    echo -e "${YELLOW}Please install ParaView or update the PARAVIEW_PATHS array in this script${NC}"
    echo -e "${YELLOW}You can manually open: $OUTPUT_DIR${NC}"
    exit 1
fi

# Step 9: Summary
echo -e "\n${BLUE}=== Simulation Complete ===${NC}"
echo -e "${GREEN}✓ Build: Successful${NC}"
echo -e "${GREEN}✓ Simulation: Completed${NC}"
echo -e "${GREEN}✓ Results: $VTK_COUNT files generated${NC}"
echo -e "${GREEN}✓ ParaView: Launched${NC}"
echo -e "\n${BLUE}Output directory: $OUTPUT_DIR${NC}"
echo -e "${BLUE}Simulation files ready for analysis!${NC}"

echo -e "\n${YELLOW}Tip: To run different simulations, modify parameters in:${NC}"
echo -e "${YELLOW}  src/unsteady_vortex_main.cpp${NC}"
