# Sweep3D Python Bindings

Python bindings for the generalized sweep library implementing "Lifted Surfacing of Generalized Sweep Volumes".

## Installation

### From Source

To build and install the Python bindings from source:

```bash
pip install .
```

Or for development:

```bash
pip install -e .
```

## Usage

### Basic Example

```python
import numpy as np
from sweep3d import GridSpec, SweepOptions, generalized_sweep

def space_time_function(point_4d):
    """
    Define a space-time implicit function.
    
    Parameters:
    -----------
    point_4d : numpy.ndarray
        4D point [x, y, z, t] where t ∈ [0, 1]
    
    Returns:
    --------
    tuple(float, numpy.ndarray)
        - Scalar value of the implicit function
        - 4D gradient [∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t]
    """
    x, y, z, t = point_4d
    
    # Example: sphere at origin with radius 0.5
    value = 0.5**2 - (x**2 + y**2 + z**2)
    gradient = np.array([-2*x, -2*y, -2*z, 0.0])
    
    return (value, gradient)

# Configure the grid
grid_spec = GridSpec()
grid_spec.resolution = [8, 8, 8]
grid_spec.bbox_min = [0.0, 0.0, 0.0]
grid_spec.bbox_max = [1.0, 1.0, 1.0]

# Configure sweep options
options = SweepOptions()
options.epsilon_env = 5e-4
options.epsilon_sil = 5e-3
options.with_adaptive_refinement = True

# Compute the sweep
result = generalized_sweep(space_time_function, grid_spec, options)

# Access the results
print(f"Envelope: {result.envelope.num_vertices()} vertices")
print(f"Sweep surface: {result.sweep_surface.num_vertices()} vertices")
```

## API Reference

### Classes

#### `GridSpec`

Specification for the 3D spatial grid used in the sweep.

**Attributes:**
- `resolution: List[int]` - Grid resolution in x, y, z directions (default: [4, 4, 4])
- `bbox_min: List[float]` - Minimum corner of bounding box (default: [-0.2, -0.2, -0.2])
- `bbox_max: List[float]` - Maximum corner of bounding box (default: [1.2, 1.2, 1.2])

#### `SweepOptions`

Options for controlling the generalized sweep operation.

**Attributes:**
- `epsilon_env: float` - Tolerance for envelope grid refinement (default: 5e-4)
- `epsilon_sil: float` - Tolerance for silhouette grid refinement (default: 5e-3)
- `max_split: int` - Maximum number of splits during refinement (default: unlimited)
- `with_insideness_check: bool` - Enable insideness checks (default: True)
- `with_snapping: bool` - Enable vertex snapping (default: True)
- `cyclic: bool` - Whether trajectory is cyclic (default: False, experimental)
- `volume_threshold: float` - Minimum volume for arrangement cells (default: 1e-5)
- `face_count_threshold: int` - Minimum face count for arrangement cells (default: 200)
- `with_adaptive_refinement: bool` - Enable adaptive grid refinement (default: True)
- `initial_time_samples: int` - Initial uniform time samples per vertex (default: 8)
- `min_tet_radius_ratio: float` - Minimum tet radius ratio (default: 1e-5)
- `min_tet_edge_length: float` - Minimum tet edge length (default: 2e-5)

#### `SweepResult`

Result of the generalized sweep operation.

**Attributes:**
- `envelope: lagrange.SurfaceMesh` - The lifted envelope mesh
- `arrangement: lagrange.SurfaceMesh` - The spatial arrangement of the envelope
- `sweep_surface: lagrange.SurfaceMesh` - The final sweep surface mesh

### Functions

#### `generalized_sweep(f, grid_spec=GridSpec(), options=SweepOptions())`

Compute the generalized sweep surface.

**Parameters:**
- `f: Callable` - Space-time implicit function taking a 4D point and returning (value, gradient)
- `grid_spec: GridSpec` - Specification of the initial spatial grid
- `options: SweepOptions` - Options controlling the sweep operation

**Returns:**
- `SweepResult` - The computed sweep result

## Space-Time Implicit Functions

A space-time implicit function defines the sweep as a level set in 4D space-time. The function should:

1. Take a 4D point `[x, y, z, t]` where `t ∈ [0, 1]` represents time
2. Return a tuple `(value, gradient)` where:
   - `value` is the scalar function value (negative outside, positive inside)
   - `gradient` is a 4D vector `[∂f/∂x, ∂f/∂y, ∂f/∂z, ∂f/∂t]`

The sweep surface is the zero level set of this function.

## Examples

See `example.py` for a complete example of computing the sweep of a moving sphere.

## License

See LICENSE.txt in the root directory.

