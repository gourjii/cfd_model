#!/usr/bin/env python3
"""
Create a ParaView-compatible time series file for the unsteady CFD simulation.
This script creates a .pvd file that ParaView can read to automatically
recognize the time series.
"""

import os
import glob
import re
from pathlib import Path

def create_pvd_file():
    # Find the output directory
    output_dir = Path("output/unsteady_vortex")
    if not output_dir.exists():
        print(f"Error: Output directory {output_dir} not found!")
        return
    
    # Find all VTK files
    vtk_files = sorted(glob.glob(str(output_dir / "unsteady_vortex_flow_flow_*.vtk")))
    
    if not vtk_files:
        print("Error: No VTK files found!")
        return
    
    print(f"Found {len(vtk_files)} VTK files")
    
    # Create PVD file content
    pvd_content = '''<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
'''
    
    for vtk_file in vtk_files:
        # Extract time step from filename
        filename = os.path.basename(vtk_file)
        match = re.search(r'_(\d{6})\.vtk$', filename)
        if match:
            time_step_int = int(match.group(1))
            # Convert to actual time (assuming dt=0.02, output every 5 steps)
            time_value = time_step_int * 0.02
            
            # Use relative path for portability
            rel_path = os.path.relpath(vtk_file, output_dir)
            
            pvd_content += f'    <DataSet timestep="{time_value:.6f}" file="{rel_path}"/>\n'
    
    pvd_content += '''  </Collection>
</VTKFile>
'''
    
    # Write PVD file
    pvd_file = output_dir / "unsteady_vortex_time_series.pvd"
    with open(pvd_file, 'w') as f:
        f.write(pvd_content)
    
    print(f"Created ParaView time series file: {pvd_file}")
    print(f"Open this file in ParaView for automatic time series recognition!")
    
    return str(pvd_file)

if __name__ == "__main__":
    create_pvd_file()
