#!/usr/bin/env python3
"""
MayaVi CFD Visualization Script
Visualizes unsteady CFD simulation results with time series animation
"""

import os
import glob
import numpy as np
from mayavi import mlab
from mayavi.sources.vtk_file_reader import VTKFileReader
import time
import argparse

class CFDVisualizer:
    def __init__(self, output_dir="output/unsteady_vortex"):
        self.output_dir = output_dir
        self.vtk_files = []
        self.current_frame = 0
        
    def find_vtk_files(self):
        """Find all VTK files in the output directory"""
        pattern = os.path.join(self.output_dir, "unsteady_vortex_flow_flow_*.vtk")
        self.vtk_files = sorted(glob.glob(pattern))
        print(f"Found {len(self.vtk_files)} VTK files")
        return len(self.vtk_files) > 0
    
    def setup_scene(self):
        """Set up the MayaVi scene"""
        # Clear any existing scene
        mlab.clf()
        
        # Set up the figure
        fig = mlab.figure(size=(1200, 800), bgcolor=(0.1, 0.1, 0.1))
        fig.name = "CFD Unsteady Vortex Flow"
        
        return fig
    
    def load_and_visualize_frame(self, frame_idx):
        """Load and visualize a specific frame"""
        if frame_idx >= len(self.vtk_files):
            return False
            
        vtk_file = self.vtk_files[frame_idx]
        print(f"Loading frame {frame_idx + 1}/{len(self.vtk_files)}: {os.path.basename(vtk_file)}")
        
        # Clear previous data
        mlab.clf()
        
        # Load VTK file
        src = mlab.pipeline.open(vtk_file)
        
        # Create velocity magnitude visualization
        velocity_mag = mlab.pipeline.extract_vector_norm(src, name='velocity_magnitude')
        
        # Create surface plot colored by velocity magnitude
        surface = mlab.pipeline.surface(velocity_mag, colormap='jet', opacity=0.8)
        surface.actor.property.interpolation = 'gouraud'
        
        # Add velocity vectors (subsampled for clarity)
        vectors = mlab.pipeline.vectors(src, scale_factor=0.1, mask_points=20)
        vectors.glyph.color_mode = 'color_by_scalar'
        vectors.glyph.glyph.scale_factor = 0.1
        
        # Add streamlines
        streamlines = mlab.pipeline.streamline(src, seedtype='plane')
        streamlines.streamline_type = 'line'
        streamlines.seed.widget.origin = [-4.0, 0.0, 0.0]
        streamlines.seed.widget.point1 = [-4.0, 2.0, 0.0]
        streamlines.seed.widget.point2 = [-4.0, -2.0, 0.0]
        streamlines.seed.widget.resolution = 10
        
        # Add obstacle outline (if we can detect it from low velocity regions)
        # This is a simple approach - in practice you might want to add the obstacle geometry separately
        
        # Set up the view
        mlab.view(azimuth=0, elevation=90, distance=15)
        
        # Add colorbar
        mlab.colorbar(surface, title='Velocity Magnitude', orientation='vertical')
        
        # Add title with time information
        time_value = frame_idx * 0.02  # Based on your time step
        mlab.title(f'CFD Unsteady Flow - Time: {time_value:.3f}s', size=0.3)
        
        # Add axes
        mlab.axes(extent=[-5, 5, -2.5, 2.5, 0, 0])
        
        return True
    
    def create_static_visualization(self, frame_idx=0):
        """Create a static visualization of a single frame"""
        self.setup_scene()
        
        if not self.find_vtk_files():
            print("No VTK files found!")
            return False
            
        success = self.load_and_visualize_frame(frame_idx)
        if success:
            print(f"Static visualization created for frame {frame_idx + 1}")
            print("Close the MayaVi window when done viewing.")
            mlab.show()
        
        return success
    
    def create_animation(self, delay=0.5, save_images=False):
        """Create an animated visualization"""
        self.setup_scene()
        
        if not self.find_vtk_files():
            print("No VTK files found!")
            return False
        
        print(f"Creating animation with {len(self.vtk_files)} frames")
        print("Press 'q' in the MayaVi window to quit the animation")
        
        # Create output directory for images if requested
        if save_images:
            img_dir = os.path.join(self.output_dir, "animation_frames")
            os.makedirs(img_dir, exist_ok=True)
            print(f"Saving animation frames to: {img_dir}")
        
        try:
            for i in range(len(self.vtk_files)):
                self.load_and_visualize_frame(i)
                
                # Save frame if requested
                if save_images:
                    img_file = os.path.join(img_dir, f"frame_{i:04d}.png")
                    mlab.savefig(img_file, size=(1200, 800))
                
                # Update display
                mlab.process_ui_events()
                
                # Check if window was closed
                if not mlab.get_engine().current_scene:
                    break
                
                # Delay between frames
                time.sleep(delay)
            
            print("Animation completed!")
            
        except KeyboardInterrupt:
            print("Animation interrupted by user")
        
        return True
    
    def create_comparison_view(self, frame_indices=[0, -1]):
        """Create a side-by-side comparison of different time steps"""
        if not self.find_vtk_files():
            print("No VTK files found!")
            return False
        
        # Ensure we have valid frame indices
        frame_indices = [max(0, min(idx if idx >= 0 else len(self.vtk_files) + idx, len(self.vtk_files) - 1)) 
                        for idx in frame_indices]
        
        print(f"Creating comparison view for frames: {frame_indices}")
        
        # Create subplots
        mlab.figure(size=(1600, 600), bgcolor=(0.1, 0.1, 0.1))
        
        for i, frame_idx in enumerate(frame_indices):
            # Create subplot
            mlab.subplot(1, len(frame_indices), i)
            
            vtk_file = self.vtk_files[frame_idx]
            src = mlab.pipeline.open(vtk_file)
            
            # Create visualization
            velocity_mag = mlab.pipeline.extract_vector_norm(src, name='velocity_magnitude')
            surface = mlab.pipeline.surface(velocity_mag, colormap='jet', opacity=0.8)
            
            # Add title
            time_value = frame_idx * 0.02
            mlab.title(f'Time: {time_value:.3f}s', size=0.2)
            
            # Set view
            mlab.view(azimuth=0, elevation=90, distance=15)
        
        mlab.show()
        return True

def main():
    parser = argparse.ArgumentParser(description='Visualize CFD simulation results with MayaVi')
    parser.add_argument('--mode', choices=['static', 'animate', 'compare'], default='static',
                       help='Visualization mode')
    parser.add_argument('--frame', type=int, default=0,
                       help='Frame index for static visualization (default: 0)')
    parser.add_argument('--delay', type=float, default=0.5,
                       help='Delay between animation frames in seconds (default: 0.5)')
    parser.add_argument('--save-frames', action='store_true',
                       help='Save animation frames as PNG images')
    parser.add_argument('--output-dir', default='output/unsteady_vortex',
                       help='Directory containing VTK files')
    
    args = parser.parse_args()
    
    # Create visualizer
    visualizer = CFDVisualizer(args.output_dir)
    
    print("=== MayaVi CFD Visualizer ===")
    print(f"Mode: {args.mode}")
    print(f"Output directory: {args.output_dir}")
    
    if args.mode == 'static':
        visualizer.create_static_visualization(args.frame)
    elif args.mode == 'animate':
        visualizer.create_animation(args.delay, args.save_frames)
    elif args.mode == 'compare':
        visualizer.create_comparison_view([0, -1])  # First and last frames
    
    print("Visualization complete!")

if __name__ == "__main__":
    main()
