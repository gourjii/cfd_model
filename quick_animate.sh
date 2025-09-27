#!/bin/bash

# Quick CFD Animation Script
# Creates animated GIF from existing unsteady simulation data
# Usage: ./quick_animate.sh
# Note: Run ./run_unsteady_flow.sh first to generate data

set -e

echo "🎬 Creating CFD Animation from existing data..."

# Activate virtual environment and run animation
source cfd_venv/bin/activate && python matplotlib_cfd_viewer.py --animate --fps 4

echo "✅ Animation complete! Check output/unsteady_vortex/cfd_animation.gif"

# Open the GIF
if [[ -f "output/unsteady_vortex/cfd_animation.gif" ]]; then
    echo "🚀 Opening animation..."
    open output/unsteady_vortex/cfd_animation.gif
fi
