#!/usr/bin/env python3
"""
Example usage of the sweep3d Python bindings.

This example demonstrates how to compute a generalized sweep surface
using a space-time implicit function.
"""

import numpy as np
from sweep3d import GridSpec, SweepOptions, generalized_sweep


def sphere_sweep_function(point_4d):
    """
    Space-time implicit function for a sphere moving along a trajectory.
    
    :param point_4d: A 4D point [x, y, z, t] where t is in [0, 1]
    :type point_4d: numpy.ndarray
    :return: Tuple containing the scalar value of the implicit function and
        a 4D gradient vector [∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t]
    :rtype: tuple(float, numpy.ndarray)
    """
    x, y, z, t = point_4d
    
    # Sphere center moves along a path (e.g., circular motion)
    center_x = 0.5 + 0.2 * np.cos(2 * np.pi * t)
    center_y = 0.5 + 0.2 * np.sin(2 * np.pi * t)
    center_z = 0.5
    
    # Sphere radius
    radius = 0.3
    
    # Distance from point to sphere center at time t
    dx = x - center_x
    dy = y - center_y
    dz = z - center_z
    
    # Implicit function value (positive inside, negative outside)
    dist_sq = dx**2 + dy**2 + dz**2
    value = radius**2 - dist_sq
    
    # Gradient computation
    # ∂f/∂x = -2(x - center_x)
    grad_x = -2 * dx
    # ∂f/∂y = -2(y - center_y)
    grad_y = -2 * dy
    # ∂f/∂z = -2(z - center_z)
    grad_z = -2 * dz
    
    # ∂f/∂t = -2(x - center_x) * (-∂center_x/∂t) - 2(y - center_y) * (-∂center_y/∂t)
    dcenter_x_dt = -0.2 * 2 * np.pi * np.sin(2 * np.pi * t)
    dcenter_y_dt = 0.2 * 2 * np.pi * np.cos(2 * np.pi * t)
    grad_t = 2 * dx * dcenter_x_dt + 2 * dy * dcenter_y_dt
    
    gradient = np.array([grad_x, grad_y, grad_z, grad_t])
    
    return (value, gradient)


def main():
    # Configure the spatial grid
    grid_spec = GridSpec()
    grid_spec.resolution = [8, 8, 8]  # 8x8x8 grid
    grid_spec.bbox_min = [0.0, 0.0, 0.0]
    grid_spec.bbox_max = [1.0, 1.0, 1.0]
    
    # Configure sweep options
    options = SweepOptions()
    options.epsilon_env = 5e-4
    options.epsilon_sil = 5e-3
    options.with_adaptive_refinement = True
    options.initial_time_samples = 8
    options.with_snapping = True
    
    # Compute the generalized sweep
    print("Computing generalized sweep...")
    result = generalized_sweep(sphere_sweep_function, grid_spec, options)
    
    print("\nSweep computation completed!")
    print(f"Envelope vertices: {result.envelope.num_vertices()}")
    print(f"Envelope faces: {result.envelope.num_facets()}")
    print(f"Arrangement vertices: {result.arrangement.num_vertices()}")
    print(f"Arrangement faces: {result.arrangement.num_facets()}")
    print(f"Sweep surface vertices: {result.sweep_surface.num_vertices()}")
    print(f"Sweep surface faces: {result.sweep_surface.num_facets()}")
    
    # Access mesh data as numpy arrays
    sweep_vertices = result.sweep_surface.get_vertices()
    sweep_facets = result.sweep_surface.get_facets()
    
    print(f"\nSweep surface vertices shape: {sweep_vertices.shape}")
    print(f"Sweep surface facets shape: {sweep_facets.shape}")
    
    # Check for attributes (e.g., "time" attribute)
    if result.sweep_surface.has_attribute("time"):
        print("Sweep surface has 'time' attribute")
    
    print(f"Available attributes: {result.sweep_surface.get_attribute_names()}")
    
    # The numpy arrays can be used with other libraries like trimesh, pyvista, etc.
    # Example with trimesh (if installed):
    # import trimesh
    # mesh = trimesh.Trimesh(vertices=sweep_vertices, faces=sweep_facets)
    # mesh.export("sweep_surface.obj")


if __name__ == "__main__":
    main()

