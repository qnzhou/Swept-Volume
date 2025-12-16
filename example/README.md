## Examples

This directory contains various examples from the paper:

* `brush_stroke_blend`: Figure 19
* `brush_stroke`: Figure 5, 6
* `fertility`: Figure 1
* `letter_L_torus`: Figure 18
* `flipping_torus`: Figure 2, 17

The following additional examples are not in the paper.

* `letter_L`
* `simple`


The easiest way to run each example is using our python binding:

```sh
# Install python binding if not done already
pip install .

# Navigate to the example directory
cd example/simple

# Run the sweep3d script with the provided yaml files
python -m sweep3d sweep.yaml config.yaml

# The output meshes will be saved in the 'output' directory
```

The output meshes are saved in [.msh format](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)
and can be viewed using [Gmsh](https://gmsh.info).
