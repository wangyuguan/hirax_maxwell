# HIRAX Maxwell

This folder contains the corner-free Chunkie HIRAX dish and its FMM3DBIE
NRCCIE tests. It does not use or load any QBX code or mesh.

Dependencies are sibling folders under `/Users/yuguan/software`:

- `fmm3dbie-hirax-dev`
- `chunkie`

Run MATLAB from this folder.

## Two-dimensional curve comparison

```matlab
run('compare_chunkie_parameter_dish.m')
```

This script only starts Chunkie. It does not construct a three-dimensional
surface, quadrature corrections, or a Maxwell system. The figure overlays
the original six-piece analytic parameter curve, the polygon
vertices passed to Chunkie, and the resulting smooth Chunkie curve. Change
the `options.*` assignments at the top of the script to experiment with the
Chunkie order, input sampling, and smoothing parameters.

## Geometry-only test

```matlab
run('test_hirax_chunkie_dish_geometry_selfconv.m')
```

This compares surface orders 4, 5, 6, and 7. It does not construct
quadrature corrections or solve a Maxwell equation.

## Maxwell NRCCIE test

```matlab
run('test_hirax_chunkie_dish_maxwell_selfconv.m')
```

This solves surface orders 4, 5, and 6 with a far electric dipole and
compares the scattered electric and magnetic fields at six far probes.
Each completed order is saved to `hirax_chunkie_dish_maxwell_selfconv.mat`.

## Close-source server run

```matlab
run('run_hirax_single_leaf.m')
```

This keeps the source and wall equator in the plane `z = 0`. The source is
placed outside the actual Chunkie-smoothed right rim with gap `R/30`. The
nominal dish radius is `R = 1`, the wavelength is `10R`, and therefore
`zk = pi/5`. Far probes are retained so the run isolates close-source
resolution from close-target evaluation. The current preflight uses a
factor-two uniform subdivision around the entire rim, one cap-collar
layer, three wall-profile panels per half, no source-side cap refinement, and
surface orders 4, 6, and 8. It creates no plots and saves each completed
order separately. Re-running the script recalculates every requested order
and overwrites files with the same names. For example, order 8 is saved as
`data/run_hirax_single_leaf_order8.mat`.

## Visualization

```matlab
run('visualize_hirax_chunkie_dish.m')
```

The figure is saved as `hirax_chunkie_dish.pdf`.

## Source

- `src/hirax_chunkie_dish_plate_surfer.m`: constructs the fixed smooth
  Chunkie planform, cap, collar, and rounded wall.
- `src/plot_surfer_patch_boundaries.m`: visualization helper.
