#!/bin/bash

# Quick CFD Animation Script
# Simple wrapper to activate virtual environment and create animation

set -e

echo "🎬 Creating CFD Animation..."

# Activate virtual environment and run animation
source cfd_venv/bin/activate && python matplotlib_cfd_viewer.py --animate --fps 4

echo "✅ Animation complete! Check output/unsteady_vortex/cfd_animation.gif"

# Open the GIF
if [[ -f "output/unsteady_vortex/cfd_animation.gif" ]]; then
    echo "🚀 Opening animation..."
    open output/unsteady_vortex/cfd_animation.gif
fi
