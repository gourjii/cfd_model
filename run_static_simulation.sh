#!/bin/bash

# CFD Static Flow Simulation Runner
# This script builds, runs the static simulation, and opens ParaView

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

# ParaView paths (common locations on macOS)
PARAVIEW_PATHS=(
    "/Applications/ParaView-5.12.0.app/Contents/MacOS/paraview"
    "/Applications/ParaView-5.11.2.app/Contents/MacOS/paraview"
    "/Applications/ParaView-5.10.1.app/Contents/MacOS/paraview"
    "/Applications/ParaView.app/Contents/MacOS/paraview"
    "/usr/local/bin/paraview"
    "/opt/homebrew/bin/paraview"
)

echo -e "${BLUE}=== CFD Static Flow Simulation Runner ===${NC}"
echo -e "${BLUE}=========================================${NC}"

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

# Navigate and build
echo -e "\n${YELLOW}Building and running static simulation...${NC}"
cd "$PROJECT_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Build
cmake .. && make || {
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
}

# Run simulation
echo -e "\n${YELLOW}Running static flow simulation...${NC}"
./bin/potential_flow_solver || {
    echo -e "${RED}✗ Simulation failed${NC}"
    exit 1
}

# Find latest results
LATEST_FLOW=$(find "$OUTPUT_DIR" -name "*potential_flow*.vtk" -type f | sort | tail -1)
LATEST_OBSTACLE=$(find "$OUTPUT_DIR" -name "*thick_obstacle*.vtk" -type f | sort | tail -1)

if [[ -n "$LATEST_FLOW" ]]; then
    echo -e "${GREEN}✓ Simulation completed successfully${NC}"
    echo -e "${GREEN}✓ Results: $(basename "$LATEST_FLOW")${NC}"
    
    # Launch ParaView
    PARAVIEW_PATH=$(find_paraview)
    if [[ $? -eq 0 ]]; then
        echo -e "\n${BLUE}Launching ParaView...${NC}"
        "$PARAVIEW_PATH" "$LATEST_FLOW" "$LATEST_OBSTACLE" &
        echo -e "${GREEN}✓ ParaView launched with flow field and obstacle${NC}"
    else
        echo -e "${YELLOW}ParaView not found. Results available at: $OUTPUT_DIR${NC}"
    fi
else
    echo -e "${RED}✗ No results found${NC}"
    exit 1
fi

echo -e "\n${BLUE}=== Static Simulation Complete ===${NC}"
