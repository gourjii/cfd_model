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

def extract_obstacle_from_vtk_data(data):
    """Extract obstacle boundary from VTK flow field data"""
    # First try to get line segments from VTK file comments
    if 'line_segments' in data and data['line_segments']:
        print(f"Found {len(data['line_segments'])} line segments from VTK comments")
        return data['line_segments']
    
    # Fallback to grid point analysis if no line segments in comments
    if data is None or 'obstacle_boundary' not in data or data['obstacle_boundary'] is None:
        return None
    
    points = data['points']
    dims = data['dimensions']
    obstacle_field = data['obstacle_boundary']
    
    # Find points where obstacle_boundary > 0.5 (marked as obstacle)
    obstacle_points = []
    
    for i in range(len(points)):
        if obstacle_field[i] > 0.5:  # This point is part of the obstacle
            obstacle_points.append([points[i][0], points[i][1]])
    
    if not obstacle_points:
        print("No obstacle points found in VTK data")
        return None
    
    obstacle_points = np.array(obstacle_points)
    print(f"Found {len(obstacle_points)} obstacle points in VTK data")
    
    # Extract line segments from obstacle points by finding the boundary structure
    line_segments = extract_line_segments_from_points(obstacle_points)
    
    print(f"Extracted {len(line_segments)} line segments from obstacle data")
    return line_segments

def extract_line_segments_from_points(obstacle_points):
    """Extract line segments from obstacle points using a generic approach"""
    if len(obstacle_points) == 0:
        return []
    
    # Generic approach: find the bounding box and create a single line
    # connecting the extreme points based on the obstacle's orientation
    
    x_coords = obstacle_points[:, 0]
    y_coords = obstacle_points[:, 1]
    
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    
    # Determine if the obstacle is more horizontal, vertical, or diagonal
    dx = x_max - x_min
    dy = y_max - y_min
    
    if dx > dy * 2:  # Mostly horizontal
        # Connect leftmost to rightmost points
        left_idx = np.argmin(x_coords)
        right_idx = np.argmax(x_coords)
        return [[obstacle_points[left_idx], obstacle_points[right_idx]]]
        
    elif dy > dx * 2:  # Mostly vertical
        # Connect topmost to bottommost points
        top_idx = np.argmax(y_coords)
        bottom_idx = np.argmin(y_coords)
        return [[obstacle_points[top_idx], obstacle_points[bottom_idx]]]
        
    else:  # Diagonal or complex shape
        # For diagonal lines, connect corner to corner
        # Find the two points that are farthest apart
        max_dist = 0
        best_pair = None
        
        for i in range(len(obstacle_points)):
            for j in range(i + 1, len(obstacle_points)):
                dist = np.sqrt((obstacle_points[i][0] - obstacle_points[j][0])**2 + 
                              (obstacle_points[i][1] - obstacle_points[j][1])**2)
                if dist > max_dist:
                    max_dist = dist
                    best_pair = (i, j)
        
        if best_pair:
            return [[obstacle_points[best_pair[0]], obstacle_points[best_pair[1]]]]
        else:
            # Fallback: connect min to max corners
            return [[[x_min, y_max], [x_max, y_min]]]

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
    
    # Read potential field data
    potential_data = None
    potential_start = None
    
    for i, line in enumerate(lines):
        if 'SCALARS potential' in line:
            # Skip the LOOKUP_TABLE line
            potential_start = i + 2
            break
    
    if potential_start is not None:
        potential_values = []
        for i in range(potential_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('SCALARS') or line.startswith('LOOKUP_TABLE'):
                break
            try:
                potential_values.append(float(line))
            except ValueError:
                break
        
        potential_data = np.array(potential_values) if potential_values else None
        if potential_data is not None:
            print(f"Read potential field with {len(potential_data)} points")
    
    # Read obstacle boundary data
    obstacle_boundary = None
    obstacle_start = None
    
    for i, line in enumerate(lines):
        if 'SCALARS obstacle_boundary' in line:
            # Skip the LOOKUP_TABLE line
            obstacle_start = i + 2
            break
    
    if obstacle_start is not None:
        obstacle_data = []
        for i in range(obstacle_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('SCALARS') or line.startswith('LOOKUP_TABLE'):
                break
            try:
                obstacle_data.append(float(line))
            except ValueError:
                break
        
        obstacle_boundary = np.array(obstacle_data) if obstacle_data else None
        if obstacle_boundary is not None:
            print(f"Read obstacle boundary field with {len(obstacle_boundary)} points")
    
    # Read pressure coefficient data (Equation 50)
    pressure_coefficient = None
    pressure_start = None
    
    for i, line in enumerate(lines):
        if 'SCALARS pressure_coefficient' in line:
            # Skip the LOOKUP_TABLE line
            pressure_start = i + 2
            break
    
    if pressure_start is not None:
        pressure_values = []
        for i in range(pressure_start, len(lines)):
            line = lines[i].strip()
            if not line or line.startswith('SCALARS') or line.startswith('LOOKUP_TABLE'):
                break
            try:
                pressure_values.append(float(line))
            except ValueError:
                break
        
        pressure_coefficient = np.array(pressure_values) if pressure_values else None
        if pressure_coefficient is not None:
            print(f"Read pressure coefficient field with {len(pressure_coefficient)} points")
    
    # Parse line segments from comments
    line_segments = []
    for line in lines:
        if line.startswith('# Segment'):
            # Parse line like: # Segment 0: (-0.5, 0.5) to (0.5, -0.5)
            match = re.search(r'# Segment \d+: \(([-\d.]+), ([-\d.]+)\) to \(([-\d.]+), ([-\d.]+)\)', line)
            if match:
                x1, y1, x2, y2 = map(float, match.groups())
                line_segments.append([[x1, y1], [x2, y2]])
    
    if line_segments:
        print(f"Read {len(line_segments)} line segments from VTK comments")
    
    return {
        'points': points,
        'dimensions': dimensions,
        'velocity': velocity_data,
        'velocity_magnitude': velocity_magnitude,
        'potential': potential_data,
        'obstacle_boundary': obstacle_boundary,
        'pressure_coefficient': pressure_coefficient,
        'line_segments': line_segments
    }

def parse_vtk_wake_vortices(filename):
    """Parse VTK POLYDATA file containing wake vortices"""
    if not os.path.exists(filename):
        return None
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find POINTS section
    points = []
    for i, line in enumerate(lines):
        if line.startswith('POINTS'):
            num_points = int(line.split()[1])
            # Read point coordinates
            for j in range(i+1, len(lines)):
                coords = lines[j].strip().split()
                if len(coords) >= 3:
                    try:
                        points.append([float(coords[0]), float(coords[1])])
                    except:
                        break
                if len(points) >= num_points:
                    break
            break
    
    return np.array(points) if points else None

def add_wake_vortices_to_plot(ax, vtk_filename):
    """Add wake vortices to existing plot"""
    wake_filename = vtk_filename.replace('_flow_flow_', '_flow_wake_')
    wake_vortices = parse_vtk_wake_vortices(wake_filename)
    if wake_vortices is not None and len(wake_vortices) > 0:
        ax.scatter(wake_vortices[:, 0], wake_vortices[:, 1], 
                  c='red', s=12, marker='o', edgecolors='darkred', 
                  linewidths=1, zorder=10, alpha=0.9)

def create_2d_visualization(data, title="CFD Visualization", save_path=None, obstacle_data=None, vtk_filename=None):
    """Create 2D visualization of CFD data with obstacle overlay - Potential + Pressure"""
    if data is None or data['points'] is None:
        return None
    
    points = data['points']
    dims = data['dimensions']
    velocity = data['velocity']
    potential = data.get('potential')
    pressure_coefficient = data.get('pressure_coefficient')
    
    # Extract 2D coordinates
    x = points[:, 0].reshape(dims[1], dims[0])  # Note: VTK uses different ordering
    y = points[:, 1].reshape(dims[1], dims[0])
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(title, fontsize=14)
    
    # Plot 1: Potential field + velocity vectors
    if potential is not None:
        potential_2d = potential.reshape(dims[1], dims[0])
        
        # Contour plot
        contour = ax1.contourf(x, y, potential_2d, levels=20, cmap='RdYlBu_r')
        ax1.set_title('Potential Field φ + Velocity Vectors')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_aspect('equal')
        plt.colorbar(contour, ax=ax1, label='Potential φ')
        
        # Add velocity vectors
        if velocity is not None:
            u = velocity[:, 0].reshape(dims[1], dims[0])
            v = velocity[:, 1].reshape(dims[1], dims[0])
            step = max(1, dims[0] // 15)
            ax1.quiver(x[::step, ::step], y[::step, ::step], 
                      u[::step, ::step], v[::step, ::step],
                      scale=None, alpha=0.7, color='black', width=0.002)
    
    # Add obstacle overlay to plot 1
    if obstacle_data is not None:
        # Plot obstacle as line segments
        for segment in obstacle_data:
            x_coords = [segment[0][0], segment[1][0]]
            y_coords = [segment[0][1], segment[1][1]]
            ax1.plot(x_coords, y_coords, 'k-', linewidth=4, label='Obstacle' if segment == obstacle_data[0] else "")
    
    # Plot 2: Pressure coefficient + velocity vectors
    if pressure_coefficient is not None:
        pressure_2d = pressure_coefficient.reshape(dims[1], dims[0])
        
        # Contour plot
        contour = ax2.contourf(x, y, pressure_2d, levels=20, cmap='coolwarm')
        ax2.set_title('Pressure Coefficient C_P + Velocity Vectors (Eq. 50)')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal')
        plt.colorbar(contour, ax=ax2, label='Pressure Coefficient C_P')
        
        # Add velocity vectors
        if velocity is not None:
            u = velocity[:, 0].reshape(dims[1], dims[0])
            v = velocity[:, 1].reshape(dims[1], dims[0])
            step = max(1, dims[0] // 15)
            ax2.quiver(x[::step, ::step], y[::step, ::step], 
                      u[::step, ::step], v[::step, ::step],
                      scale=None, alpha=0.7, color='black', width=0.002)
    
    # Add obstacle overlay to plot 2
    if obstacle_data is not None:
        # Plot obstacle as line segments
        for segment in obstacle_data:
            x_coords = [segment[0][0], segment[1][0]]
            y_coords = [segment[0][1], segment[1][1]]
            ax2.plot(x_coords, y_coords, 'k-', linewidth=4, label='Obstacle' if segment == obstacle_data[0] else "")
    
    # Add wake vortices to both plots
    if vtk_filename is not None:
        add_wake_vortices_to_plot(ax1, vtk_filename)
        add_wake_vortices_to_plot(ax2, vtk_filename)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {save_path}")
    
    return fig

def create_animation(vtk_files, output_path="cfd_animation.gif", fps=2, obstacle_data=None):
    """Create animated GIF from VTK files with side-by-side potential and pressure plots"""
    print(f"Creating animation from {len(vtk_files)} frames...")
    
    # Read first frame to set up plot
    first_data = parse_vtk_structured_grid(vtk_files[0])
    if first_data is None:
        return False
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    
    # Store colorbars to update them
    cbar1 = None
    cbar2 = None
    
    def animate(frame_idx):
        nonlocal cbar1, cbar2
        
        # Clear both axes
        ax1.clear()
        ax2.clear()
        
        # Read data for this frame
        data = parse_vtk_structured_grid(vtk_files[frame_idx])
        if data is None:
            return
        
        points = data['points']
        dims = data['dimensions']
        velocity = data['velocity']
        potential = data.get('potential')
        pressure_coefficient = data.get('pressure_coefficient')
        
        # Extract 2D coordinates
        x = points[:, 0].reshape(dims[1], dims[0])
        y = points[:, 1].reshape(dims[1], dims[0])
        
        # LEFT PLOT: Potential field + velocity vectors
        if potential is not None:
            potential_2d = potential.reshape(dims[1], dims[0])
            
            # Contour plot
            contour1 = ax1.contourf(x, y, potential_2d, levels=20, cmap='RdYlBu_r')
            if cbar1 is None:
                cbar1 = plt.colorbar(contour1, ax=ax1, label='Potential φ')
            
            # Add velocity vectors
            if velocity is not None:
                u = velocity[:, 0].reshape(dims[1], dims[0])
                v = velocity[:, 1].reshape(dims[1], dims[0])
                
                # Subsample for vector plot
                step = max(1, dims[0] // 15)
                ax1.quiver(x[::step, ::step], y[::step, ::step], 
                          u[::step, ::step], v[::step, ::step],
                          scale=None, alpha=0.7, color='black', width=0.002)
        
        ax1.set_title('Potential Field φ + Velocity Vectors')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_aspect('equal')
        
        # RIGHT PLOT: Pressure coefficient + velocity vectors
        if pressure_coefficient is not None:
            pressure_2d = pressure_coefficient.reshape(dims[1], dims[0])
            
            # Contour plot
            contour2 = ax2.contourf(x, y, pressure_2d, levels=20, cmap='coolwarm')
            if cbar2 is None:
                cbar2 = plt.colorbar(contour2, ax=ax2, label='C_P')
            
            # Add velocity vectors
            if velocity is not None:
                u = velocity[:, 0].reshape(dims[1], dims[0])
                v = velocity[:, 1].reshape(dims[1], dims[0])
                
                # Subsample for vector plot
                step = max(1, dims[0] // 15)
                ax2.quiver(x[::step, ::step], y[::step, ::step], 
                          u[::step, ::step], v[::step, ::step],
                          scale=None, alpha=0.7, color='black', width=0.002)
        
        ax2.set_title('Pressure Coefficient C_P + Velocity (Eq. 50)')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_aspect('equal')
        
        # Add obstacle overlay to both plots
        if obstacle_data is not None:
            for segment in obstacle_data:
                x_coords = [segment[0][0], segment[1][0]]
                y_coords = [segment[0][1], segment[1][1]]
                ax1.plot(x_coords, y_coords, 'k-', linewidth=4)
                ax2.plot(x_coords, y_coords, 'k-', linewidth=4)
        
        # Add wake vortices overlay to both plots
        wake_filename = vtk_files[frame_idx].replace('_flow_flow_', '_flow_wake_')
        wake_vortices = parse_vtk_wake_vortices(wake_filename)
        if wake_vortices is not None and len(wake_vortices) > 0:
            ax1.scatter(wake_vortices[:, 0], wake_vortices[:, 1], 
                       c='red', s=12, marker='o', edgecolors='darkred', 
                       linewidths=1, zorder=10, alpha=0.9, label='Wake Vortices')
            ax2.scatter(wake_vortices[:, 0], wake_vortices[:, 1], 
                       c='red', s=12, marker='o', edgecolors='darkred', 
                       linewidths=1, zorder=10, alpha=0.9)
        
        time_value = frame_idx * 0.02
        fig.suptitle(f'CFD Unsteady Flow - Time: {time_value:.3f}s', fontsize=16, fontweight='bold')
    
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
    
    # Find VTK files (try different patterns)
    patterns = [
        os.path.join(args.output_dir, "unsteady_vortex_flow_flow_*.vtk"),  # Unsteady pattern
        os.path.join(args.output_dir, "*potential_flow*.vtk"),             # Static pattern
        os.path.join(args.output_dir, "*.vtk")                             # General pattern
    ]
    
    vtk_files = []
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            vtk_files = sorted(files)
            break
    
    if not vtk_files:
        print(f"No VTK files found in {args.output_dir}")
        return
    
    print(f"Found {len(vtk_files)} VTK files")
    
    # Extract obstacle data from the first VTK file (all should have the same obstacle)
    obstacle_data = None
    if vtk_files:
        first_vtk_data = parse_vtk_structured_grid(vtk_files[0])
        obstacle_data = extract_obstacle_from_vtk_data(first_vtk_data)
    
    if obstacle_data is None:
        print("No obstacle boundary found in VTK data - visualizing flow field only")
    
    if args.animate:
        # Create animation
        output_path = os.path.join(args.output_dir, "cfd_animation.gif")
        create_animation(vtk_files, output_path, args.fps, obstacle_data)
    elif args.save_frames:
        # Save all frames as images
        img_dir = os.path.join(args.output_dir, "frame_images")
        os.makedirs(img_dir, exist_ok=True)
        
        for i, vtk_file in enumerate(vtk_files):
            data = parse_vtk_structured_grid(vtk_file)
            time_value = i * 0.02
            title = f"CFD Flow - Time: {time_value:.3f}s"
            
            img_path = os.path.join(img_dir, f"frame_{i:04d}.png")
            fig = create_2d_visualization(data, title, img_path, obstacle_data, vtk_file)
            if fig:
                plt.close(fig)
    else:
        # Single frame visualization
        frame_idx = max(0, min(args.frame, len(vtk_files) - 1))
        vtk_file = vtk_files[frame_idx]
        
        data = parse_vtk_structured_grid(vtk_file)
        time_value = frame_idx * 0.02
        title = f"CFD Flow - Time: {time_value:.3f}s (Frame {frame_idx+1}/{len(vtk_files)})"
        
        fig = create_2d_visualization(data, title, None, obstacle_data, vtk_file)
        if fig:
            plt.show()

if __name__ == "__main__":
    main()
