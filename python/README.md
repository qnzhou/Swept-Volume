# Sweep3D Python Bindings

Python bindings for the generalized sweep library implementing "Lifted Surfacing of Generalized Sweep Volumes".

## Installation

```bash
pip install .
# Or for development
pip install -e .
```

## Quick Start

```python
import numpy as np
from sweep3d import GridSpec, SweepOptions, generalized_sweep

def space_time_function(point_4d):
    """Space-time implicit function.
    
    :param point_4d: 4D point [x, y, z, t] where t ∈ [0, 1]
    :type point_4d: numpy.ndarray
    :return: Tuple of (value, gradient) where gradient is [∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t]
    :rtype: tuple(float, numpy.ndarray)
    """
    x, y, z, t = point_4d
    value = 0.5**2 - (x**2 + y**2 + z**2)
    gradient = np.array([-2*x, -2*y, -2*z, 0.0])
    return (value, gradient)

# Configure and compute
grid_spec = GridSpec()
grid_spec.resolution = [8, 8, 8]
grid_spec.bbox_min = [0.0, 0.0, 0.0]
grid_spec.bbox_max = [1.0, 1.0, 1.0]

options = SweepOptions()
options.epsilon_env = 5e-4
options.epsilon_sil = 5e-3

result = generalized_sweep(space_time_function, grid_spec, options)
print(f"Sweep surface: {result.sweep_surface.num_vertices()} vertices")
```

## API Reference

### `generalized_sweep(f, grid_spec, options)`

Compute the generalized sweep surface.

**Parameters:**
- `f` (Callable): Space-time implicit function taking a 4D point [x, y, z, t] and returning (value, gradient)
- `grid_spec` (GridSpec, optional): Specification of the initial spatial grid
- `options` (SweepOptions, optional): Options controlling the sweep operation

**Returns:**
- `SweepResult`: Result containing envelope, arrangement, and sweep_surface meshes

### `GridSpec`

**Attributes:**
- `resolution` (List[int]): Grid resolution [nx, ny, nz] (default: [4, 4, 4])
- `bbox_min` (List[float]): Bounding box minimum (default: [-0.2, -0.2, -0.2])
- `bbox_max` (List[float]): Bounding box maximum (default: [1.2, 1.2, 1.2])

### `SweepOptions`

**Key Attributes:**
- `epsilon_env` (float): Envelope tolerance (default: 5e-4)
- `epsilon_sil` (float): Silhouette tolerance (default: 5e-3)
- `with_adaptive_refinement` (bool): Enable adaptive refinement (default: True)
- `with_snapping` (bool): Enable vertex snapping (default: True)
- `initial_time_samples` (int): Initial time samples per vertex (default: 8)
- `cyclic` (bool): Cyclic trajectory (default: False, experimental)

### `SweepResult`

**Attributes:**
- `envelope` (lagrange.SurfaceMesh): Lifted envelope mesh
- `arrangement` (lagrange.SurfaceMesh): Spatial arrangement mesh
- `sweep_surface` (lagrange.SurfaceMesh): Final sweep surface mesh

## Space-Time Functions

Define sweeps as 4D implicit functions:
- Input: 4D point [x, y, z, t] where t ∈ [0, 1]
- Output: (value, gradient) where gradient = [∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t]
- Surface: Zero level set of the function


