# Lifted Surfacing of Generalized Sweep Volumes

This code implements the ACM SIGGRAPH ASIA 2025 paper: Lifted Surfacing of Generalized Sweep Volumes

<img width="1600" alt="ball-rolling" src="https://github.com/user-attachments/assets/4751b437-091b-4626-976e-3c0a74132838" />

>Top: A wire-like ball rolls forward while offsetting, changing its genus from 41 to 29. Bottom: The time-colored sweep boundary and sharp creases (2nd row), in transparency (3rd row, creases hidden for clarity), and a cut-away view (bottom row).

Given any sweep represented as a smooth time-varying implicit function satisfying a genericity assumption, this algorithm produces a watertight and intersection-free surface that faithfully approximates the geometric and topological features.

## Build

### C++ Build

Use the following command to build: 

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
The program `general_sweep` will be generated in the build file. 

### Python Bindings

Python bindings are available for easy integration with Python workflows:

```bash
# Install the Python package
pip install .

# Or install in development mode
pip install -e .

# Quick build script
./build_python.sh
```

See [`python/README.md`](python/README.md) for detailed Python API documentation and examples.

### Dependency

Currently, all third-party libraries are now open-sourced.

## Usage

### Python API

For Python users, the library provides a simple API:

```python
import numpy as np
from sweep3d import GridSpec, SweepOptions, generalized_sweep

def my_space_time_function(point_4d):
    """Define a space-time implicit function."""
    x, y, z, t = point_4d
    # Example: static sphere
    value = 0.5**2 - (x**2 + y**2 + z**2)
    gradient = np.array([-2*x, -2*y, -2*z, 0.0])
    return (value, gradient)

# Configure and compute
grid_spec = GridSpec()
grid_spec.resolution = [8, 8, 8]
options = SweepOptions()
result = generalized_sweep(my_space_time_function, grid_spec, options)

# Access results
print(f"Sweep surface: {result.sweep_surface.num_vertices()} vertices")
```

See [`python/example.py`](python/example.py) for a complete working example and [`python/README.md`](python/README.md) for detailed API documentation.

### Command Line Tool

To use the `general_sweep` tool, you must provide an initial grid file and output path as required arguments, along with any desired options. The trajectory functions are defined in `trajectory.h`. 

```bash
./general_sweep <grid> <output> [OPTIONS]
```


### 4D Implicit Function Framework: 

The input of this program is any generalized sweep that is represented by a smooth 4D implicit function. Currently, we provide a series of pre-defined functions which include all paper examples. You can specify an implicit function file or use one of the predefined function names. Unfortunately, we do not provide a GUI or a user-friendly tool for defining such functions at this moment. If you want to specify your own sweep, please refer to this [repository](https://github.com/adobe-research/space-time-functions) and specific use cases in `trajectory.h` for details.

### Positional Arguments

- `grid` : The path to the initial grid file that will be used for grid generation. This file can be either a `.msh` or `.json` file. When a `.json` file is provided, it will be converted to a `.msh` file internally. The format to specify an initial grid can be found [here](https://github.com/Jurwen/Swept-Volume/blob/main/data/test/grid_1.json). For this file format, you only need to specify the bounding box of the initial grid and how "dense" this initial grid is, which is controlled by the `resolution` parameter. The higher the resolution, the denser the initial mesh.
- `output` : The output directory path where all generated files will be saved. The tool will create this directory if it doesn't exist. `.obj` files only contain the mesh, while `.msh` files have additional time info. Output files include:
  - `envelope.msh`, `envelope.obj` : The mesh before arrangement (if `SAVE_CONTOUR` is enabled)
  - (`0.obj`, `1.obj`, ...), (`0.msh`, `1.msh`, ...) : Separated cell components with 0 winding number.
  - `mesh.obj`, `mesh.msh` : Combined cell componnets with 0 winding number.
  - `features.json` : Feature lines (first slot) and points (second slot) 

### Options

- `-h, --help` : Show the help message and exit the program.
- `-f, --function <file>` : Specify an implicit function file or predefined function name. Can be:
  - A predefined function name (e.g., `fertility_v4`, `kitten_dog`, `letter_L_blend`, `ball_genus_roll`, `tangle_chair_S`, `star_S`, etc.)
  - See `trajectory.h` for all available predefined functions
- `--ee, --epsilon-env <value>` : Set the envelope threshold (default: 0.0005). Lower values produce finer grids for the sweep function. This is a DOUBLE value that controls the envelope's precision. This parameter corresponds to the `epsilon_env` variable stated in the paper
- `--es, --epsilon-sil <value>` : Set the silhouette threshold (default: 0.005). Lower values produce finer grids for the silhouette function (sweep function taking partial derivative in time).This is a DOUBLE value that controls trajectory's' precision. This parameter corresponds to the `epsilon_sil` variable stated in the paper
- `-i, --inside-check`: Whether the grid generation will do a full grid refinement considering regions of the envelope that is inside the sweep. Turning this on will generate a complete smooth envelope with fine details for inside regions. Note that this tag may add considerable amount of time and memory. 

## Example:

The following is an example of how to use the `general_sweep` tool with common options:

```bash
./general_sweep ../data/test/grid_1.json ../output/brush_stroke_example --ee 0.0005 --es 0.005 -f brush_stroke_blending
```

### Parameter Breakdown:

This example command demonstrates how to generate a swept volume using the `brush_stroke_blending` function. Let's examine each parameter:

#### Required Positional Arguments:
- **`../data/test/grid_1.json`** : The initial grid file that defines the starting tetrahedral mesh structure. This JSON file contains the initial spatial discretization that will be refined during the swept volume computation.

- **`../output/brush_stroke_example`** : The output directory where all generated files will be saved. The tool will create this directory if it doesn't exist.

#### Optional Parameters:
- **`--ee 0.0005`** : Sets the **envelope threshold** to 0.0005. This parameter controls how finely the algorithm subdivides the initial grid based on the implicit function's complexity. A small value (like 0.0005) means:
  - High precision in capturing surface details
  - More computational time and memory usage
  - Better preservation of sharp features and fine geometric details
  - The algorithm will subdivide grid cells more aggressively where the function changes rapidly

- **`--es 0.005`** : Sets the **silhouette threshold** to 0.005. This parameter controls the precision of trajectory processing, which is crucial for:
  - Temporal discretization of the 4D sweep
  - Determining how finely to sample the 4D grid based on the complexity of the trajectory
  - Capturing temporal variations in the implicit function

- **`-f brush_stroke_blending`** : Specifies the **implicit function** to use. The `brush_stroke_blending` function represents a particular sweep pattern that:
  - Creates a brush stroke-like swept volume with blending effects
  - Demonstrates complex topological changes during the sweep
  - Shows how multiple geometric elements can be blended together in the swept volume, as the base brush morphs between one sphere and two spheres.
  - Is one of the predefined functions available in `trajectory.h`

### Expected Output:
When this command runs successfully, you'll find in the `../output/brush_stroke_example/` directory:
- **`0.obj`, `1.obj`, ...** : Individual mesh components representing different parts of the swept envelope with 0 winding number
- **`0.msh`, `1.msh`, ...** : Containing per-vertex info of `time`
- **`mesh.obj`** : Combined mesh components of the swept envelope with 0 winding number
- **`mesh.msh`** : Containing per-vertex info of `time`
- **`envelope.obj`** : The intermediate envelope mesh before pruning
- **`envelope.msh`** : Containing per-vertex info of `time` and `is_regular`
- **`features.json`** : Feature lines and points that capture the topological structure of the swept volume

This example showcases how different parameters should work in an actual use case.
