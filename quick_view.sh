#!/bin/bash

# Quick CFD Viewer Script
# Simple wrapper to view a specific frame

set -e

FRAME=${1:-5}  # Default to frame 5 if no argument provided

echo "👁️  Viewing CFD Frame $FRAME..."

# Activate virtual environment and show frame
source cfd_venv/bin/activate && python matplotlib_cfd_viewer.py --frame $FRAME

echo "✅ Viewer complete!"
