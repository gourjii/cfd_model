#!/usr/bin/env python3
"""
Simple MayaVi CFD Viewer
A robust viewer for CFD VTK files with error handling
"""

import os
import glob
import numpy as np
from mayavi import mlab
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import argparse

def read_vtk_structured_grid(filename):
    """Read VTK structured grid file and extract data"""
    print(f"Reading: {filename}")
    
    # Read VTK file
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    
    # Get the data
    output = reader.GetOutput()
    
    # Get points
    points = output.GetPoints()
    if points is None:
        print("No points found in VTK file")
        return None
        
    points_array = vtk_to_numpy(points.GetData())
    
    # Get dimensions
    dims = output.GetDimensions()
    print(f"Grid dimensions: {dims}")
    
    # Get point data
    point_data = output.GetPointData()
    
    # Extract velocity vectors
    velocity = None
    velocity_magnitude = None
    
    for i in range(point_data.GetNumberOfArrays()):
        array_name = point_data.GetArrayName(i)
        array = point_data.GetArray(i)
        print(f"Found array: {array_name} with {array.GetNumberOfComponents()} components")
        
        if array_name == 'velocity' and array.GetNumberOfComponents() == 3:
            velocity = vtk_to_numpy(array)
        elif array_name == 'velocity_magnitude':
            velocity_magnitude = vtk_to_numpy(array)
    
    # If we have velocity but no magnitude, calculate it
    if velocity is not None and velocity_magnitude is None:
        velocity_magnitude = np.sqrt(np.sum(velocity**2, axis=1))
    
    return {
        'points': points_array,
        'dimensions': dims,
        'velocity': velocity,
        'velocity_magnitude': velocity_magnitude
    }

def visualize_cfd_data(data, title="CFD Visualization"):
    """Create visualization from CFD data"""
    if data is None:
        return False
    
    points = data['points']
    dims = data['dimensions']
    velocity_mag = data['velocity_magnitude']
    velocity = data['velocity']
    
    # Reshape data to grid
    x = points[:, 0].reshape(dims, order='F')
    y = points[:, 1].reshape(dims, order='F')
    z = points[:, 2].reshape(dims, order='F')
    
    if velocity_mag is not None:
        scalar_data = velocity_mag.reshape(dims, order='F')
    else:
        print("No velocity magnitude data found")
        return False
    
    # Clear scene
    mlab.clf()
    
    # Create surface plot
    surf = mlab.mesh(x, y, z, scalars=scalar_data, colormap='jet')
    
    # Add velocity vectors if available (subsampled)
    if velocity is not None:
        # Subsample for clarity
        step = max(1, dims[0] // 20)  # Show every 20th point
        indices = np.arange(0, len(points), step)
        
        if len(indices) > 0:
            pos = points[indices]
            vel = velocity[indices]
            
            # Scale vectors for visibility
            scale = 0.1
            mlab.quiver3d(pos[:, 0], pos[:, 1], pos[:, 2],
                         vel[:, 0], vel[:, 1], vel[:, 2],
                         scale_factor=scale, color=(1, 1, 1))
    
    # Add colorbar
    mlab.colorbar(surf, title='Velocity Magnitude', orientation='vertical')
    
    # Set view
    mlab.view(azimuth=0, elevation=90, distance='auto')
    
    # Add title
    mlab.title(title, size=0.3)
    
    # Add outline
    mlab.outline()
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Simple MayaVi CFD Viewer')
    parser.add_argument('--output-dir', default='output/unsteady_vortex',
                       help='Directory containing VTK files')
    parser.add_argument('--frame', type=int, default=0,
                       help='Frame index to visualize (default: 0)')
    parser.add_argument('--animate', action='store_true',
                       help='Animate through all frames')
    parser.add_argument('--delay', type=float, default=1.0,
                       help='Delay between frames in animation (seconds)')
    
    args = parser.parse_args()
    
    # Find VTK files
    pattern = os.path.join(args.output_dir, "unsteady_vortex_flow_flow_*.vtk")
    vtk_files = sorted(glob.glob(pattern))
    
    if not vtk_files:
        print(f"No VTK files found in {args.output_dir}")
        return
    
    print(f"Found {len(vtk_files)} VTK files")
    
    # Set up MayaVi
    mlab.figure(size=(1200, 800), bgcolor=(0.1, 0.1, 0.1))
    
    if args.animate:
        print("Starting animation...")
        print("Close the MayaVi window to stop")
        
        try:
            import time
            for i, vtk_file in enumerate(vtk_files):
                data = read_vtk_structured_grid(vtk_file)
                time_value = i * 0.02  # Based on your time step
                title = f"CFD Flow - Time: {time_value:.3f}s (Frame {i+1}/{len(vtk_files)})"
                
                if visualize_cfd_data(data, title):
                    mlab.process_ui_events()
                    time.sleep(args.delay)
                else:
                    print(f"Failed to visualize frame {i}")
                    break
                    
        except KeyboardInterrupt:
            print("Animation stopped by user")
    else:
        # Single frame
        frame_idx = max(0, min(args.frame, len(vtk_files) - 1))
        vtk_file = vtk_files[frame_idx]
        
        data = read_vtk_structured_grid(vtk_file)
        time_value = frame_idx * 0.02
        title = f"CFD Flow - Time: {time_value:.3f}s (Frame {frame_idx+1}/{len(vtk_files)})"
        
        if visualize_cfd_data(data, title):
            print("Visualization created. Close window when done.")
            mlab.show()
        else:
            print("Failed to create visualization")

if __name__ == "__main__":
    main()
