#!/usr/bin/env python3
"""
Matplotlib CFD Viewer
A reliable viewer for CFD VTK files using matplotlib
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import argparse
import re

def parse_vtk_structured_grid(filename):
    """Parse VTK structured grid file manually"""
    print(f"Reading: {filename}")
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    dimensions = None
    points_start = None
    point_data_start = None
    
    for i, line in enumerate(lines):
        if line.startswith('DIMENSIONS'):
            dims = [int(x) for x in line.split()[1:]]
            dimensions = dims
            print(f"Dimensions: {dimensions}")
        elif line.startswith('POINTS'):
            points_start = i + 1
            num_points = int(line.split()[1])
            print(f"Number of points: {num_points}")
        elif line.startswith('POINT_DATA'):
            point_data_start = i + 1
            break
    
    if dimensions is None or points_start is None:
        print("Could not parse VTK file structure")
        return None
    
    # Read points
    points = []
    for i in range(points_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('FIELD') or line.startswith('POINT_DATA') or line.startswith('VECTORS'):
            break
        coords = [float(x) for x in line.split()]
        for j in range(0, len(coords), 3):
            if j + 2 < len(coords):
                points.append([coords[j], coords[j+1], coords[j+2]])
    
    points = np.array(points)
    print(f"Read {len(points)} points")
    
    # Read velocity data
    velocity_data = []
    velocity_start = None
    
    for i, line in enumerate(lines):
        if 'VECTORS velocity' in line:
            velocity_start = i + 1
            break
    
    if velocity_start is not None:
        for i in range(velocity_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('SCALARS') or line.startswith('LOOKUP_TABLE'):
                break
            velocities = [float(x) for x in line.split()]
            for j in range(0, len(velocities), 3):
                if j + 2 < len(velocities):
                    velocity_data.append([velocities[j], velocities[j+1], velocities[j+2]])
    
    velocity_data = np.array(velocity_data) if velocity_data else None
    
    if velocity_data is not None:
        print(f"Read {len(velocity_data)} velocity vectors")
        velocity_magnitude = np.sqrt(np.sum(velocity_data**2, axis=1))
    else:
        print("No velocity data found")
        velocity_magnitude = None
    
    return {
        'points': points,
        'dimensions': dimensions,
        'velocity': velocity_data,
        'velocity_magnitude': velocity_magnitude
    }

def create_2d_visualization(data, title="CFD Visualization", save_path=None):
    """Create 2D visualization of CFD data"""
    if data is None or data['points'] is None:
        return None
    
    points = data['points']
    dims = data['dimensions']
    velocity = data['velocity']
    velocity_mag = data['velocity_magnitude']
    
    # Extract 2D coordinates
    x = points[:, 0].reshape(dims[1], dims[0])  # Note: VTK uses different ordering
    y = points[:, 1].reshape(dims[1], dims[0])
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(title, fontsize=14)
    
    # Plot 1: Velocity magnitude contour
    if velocity_mag is not None:
        vel_mag_2d = velocity_mag.reshape(dims[1], dims[0])
        
        # Contour plot
        contour = ax1.contourf(x, y, vel_mag_2d, levels=20, cmap='jet')
        ax1.set_title('Velocity Magnitude')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_aspect('equal')
        plt.colorbar(contour, ax=ax1, label='Velocity Magnitude')
        
        # Add contour lines
        ax1.contour(x, y, vel_mag_2d, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    
    # Plot 2: Velocity vectors
    if velocity is not None:
        u = velocity[:, 0].reshape(dims[1], dims[0])
        v = velocity[:, 1].reshape(dims[1], dims[0])
        
        # Subsample for vector plot
        step = max(1, dims[0] // 20)
        x_sub = x[::step, ::step]
        y_sub = y[::step, ::step]
        u_sub = u[::step, ::step]
        v_sub = v[::step, ::step]
        
        # Vector plot with velocity magnitude as background
        if velocity_mag is not None:
            vel_mag_2d = velocity_mag.reshape(dims[1], dims[0])
            im = ax2.imshow(vel_mag_2d, extent=[x.min(), x.max(), y.min(), y.max()], 
                           origin='lower', cmap='jet', alpha=0.7)
            plt.colorbar(im, ax=ax2, label='Velocity Magnitude')
        
        # Add velocity vectors
        ax2.quiver(x_sub, y_sub, u_sub, v_sub, scale=None, alpha=0.8, color='white', width=0.003)
        ax2.set_title('Velocity Vectors')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig

def create_animation(vtk_files, output_path="cfd_animation.gif", fps=2):
    """Create animated GIF from VTK files"""
    print(f"Creating animation from {len(vtk_files)} frames...")
    
    # Read first frame to set up plot
    first_data = parse_vtk_structured_grid(vtk_files[0])
    if first_data is None:
        return False
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    def animate(frame_idx):
        ax.clear()
        
        # Read data for this frame
        data = parse_vtk_structured_grid(vtk_files[frame_idx])
        if data is None:
            return
        
        points = data['points']
        dims = data['dimensions']
        velocity_mag = data['velocity_magnitude']
        velocity = data['velocity']
        
        # Extract 2D coordinates
        x = points[:, 0].reshape(dims[1], dims[0])
        y = points[:, 1].reshape(dims[1], dims[0])
        
        if velocity_mag is not None:
            vel_mag_2d = velocity_mag.reshape(dims[1], dims[0])
            
            # Contour plot
            contour = ax.contourf(x, y, vel_mag_2d, levels=20, cmap='jet')
            
            # Add velocity vectors
            if velocity is not None:
                u = velocity[:, 0].reshape(dims[1], dims[0])
                v = velocity[:, 1].reshape(dims[1], dims[0])
                
                # Subsample for vector plot
                step = max(1, dims[0] // 15)
                x_sub = x[::step, ::step]
                y_sub = y[::step, ::step]
                u_sub = u[::step, ::step]
                v_sub = v[::step, ::step]
                
                ax.quiver(x_sub, y_sub, u_sub, v_sub, scale=None, alpha=0.8, 
                         color='white', width=0.002)
        
        time_value = frame_idx * 0.02
        ax.set_title(f'CFD Unsteady Flow - Time: {time_value:.3f}s', fontsize=14)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(vtk_files), 
                                  interval=1000//fps, repeat=True)
    
    # Save as GIF
    print(f"Saving animation to: {output_path}")
    anim.save(output_path, writer='pillow', fps=fps)
    print("Animation saved!")
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Matplotlib CFD Viewer')
    parser.add_argument('--output-dir', default='output/unsteady_vortex',
                       help='Directory containing VTK files')
    parser.add_argument('--frame', type=int, default=0,
                       help='Frame index to visualize (default: 0)')
    parser.add_argument('--animate', action='store_true',
                       help='Create animated GIF')
    parser.add_argument('--save-frames', action='store_true',
                       help='Save individual frame images')
    parser.add_argument('--fps', type=int, default=2,
                       help='Frames per second for animation (default: 2)')
    
    args = parser.parse_args()
    
    # Find VTK files
    pattern = os.path.join(args.output_dir, "unsteady_vortex_flow_flow_*.vtk")
    vtk_files = sorted(glob.glob(pattern))
    
    if not vtk_files:
        print(f"No VTK files found in {args.output_dir}")
        return
    
    print(f"Found {len(vtk_files)} VTK files")
    
    if args.animate:
        # Create animation
        output_path = os.path.join(args.output_dir, "cfd_animation.gif")
        create_animation(vtk_files, output_path, args.fps)
    elif args.save_frames:
        # Save all frames as images
        img_dir = os.path.join(args.output_dir, "frame_images")
        os.makedirs(img_dir, exist_ok=True)
        
        for i, vtk_file in enumerate(vtk_files):
            data = parse_vtk_structured_grid(vtk_file)
            time_value = i * 0.02
            title = f"CFD Flow - Time: {time_value:.3f}s"
            
            img_path = os.path.join(img_dir, f"frame_{i:04d}.png")
            fig = create_2d_visualization(data, title, img_path)
            if fig:
                plt.close(fig)
    else:
        # Single frame visualization
        frame_idx = max(0, min(args.frame, len(vtk_files) - 1))
        vtk_file = vtk_files[frame_idx]
        
        data = parse_vtk_structured_grid(vtk_file)
        time_value = frame_idx * 0.02
        title = f"CFD Flow - Time: {time_value:.3f}s (Frame {frame_idx+1}/{len(vtk_files)})"
        
        fig = create_2d_visualization(data, title)
        if fig:
            plt.show()

if __name__ == "__main__":
    main()
