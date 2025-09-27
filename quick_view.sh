#!/bin/bash

# Quick CFD Frame Viewer Script
# Views a specific frame from unsteady simulation data
# Usage: ./quick_view.sh [frame_number]
# Example: ./quick_view.sh 3
# Note: Run ./run_unsteady_flow.sh first to generate data

set -e

FRAME=${1:-5}  # Default to frame 5 if no argument provided

echo "👁️  Viewing CFD Frame $FRAME from existing data..."

# Activate virtual environment and show frame
source cfd_venv/bin/activate && python matplotlib_cfd_viewer.py --frame $FRAME

echo "✅ Viewer complete!"
