# Dish geometry

Imported from the shared `Dropbox/hirax/geometry/dish` directory on 2026-10-12.
The CAD generator and numerical meshing/smoothing algorithms are retained;
dependency and output paths have been adapted to this repository.

The reflector is a closed thin shell with a 6 m aperture and a 1 cm wall.
Its inner surface is the paraboloid `z = r^2/(4*f) - f`, with focus at the
origin. The outer surface approximates a normal offset, and a narrow rim
joins the two surfaces. There is no disk across the aperture. This thickness
was added to the original zero-thickness CST reflector model.

## Files

- `cad/`: Python generator and four CAD variants, each with `.geo`, `.step`,
  and analytic `.json` metadata. `fd021` and `fd023` mean f/D = 0.21 and 0.23.
  `_m` uses metres; `_t1` uses centimetres for the same physical dimensions.
- `mesh/run_dish_mesh.m`: STEP meshing, geometric diagnostics, and fmm3dbie
  `.go3` / `.gidmsh` export. Generated files go into `mesh/out/` and are ignored.
- `preview_dish.m`: load a `.go3` file through fmm3dbie and plot the surface
  and high-order patch boundaries.
- `experimental/`: unfinished smoothing workflows; see its README before use.

This module currently contains geometry tools, not a dish scattering solver.
The imported meshes have not been validated for a converged scattering solve.

## View the existing mesh

From the repository root in MATLAB:

```matlab
addpath('src/dish');
setup_dish;
[S,fig] = preview_dish();
% Or supply any exported mesh explicitly:
% [S,fig] = preview_dish('/path/to/dish_order4.go3');
```

The default mesh is
`data/dish/dish_fd023_m_rimfree/dish_fd023_m_rimfree_order4.go3`.
A local copy of this mesh, its scaffold, and its report was imported from the
shared folder. It contains 5,966 order-4 patches and 89,490 nodes. Generated
mesh data are excluded from Git. On another machine, copy that directory from
the shared folder's `mesh/out/`, generate a mesh as below, or pass a mesh path
to `preview_dish`. `DISH_MESH_ROOT` can override `data/dish`.
Preview `.fig` and `.png` files go into the ignored `results/` directory.

The adaptive f/D = 0.23 candidate produced on 2026-10-13 is tracked with
Git LFS at
`data/dish/dish_fd023_m_adaptive_ad_as0_sig050w_cad12_p04_eps1e03/`.
It contains 280,296 order-4 patches and 4,204,440 nodes. Its surface-validity
check passed, but adaptive refinement stopped at depth 3 with 12 patches still
above the requested tolerance, so it remains a candidate rather than a
converged reference mesh. Install Git LFS before cloning or run `git lfs pull`
to materialize the 1.18 GB `.go3` file.

`setup_dish` uses fmm3dbie already on the MATLAB path, or looks for
`fmm3dbie-hirax-dev` beside this repository. Set `NRCCIE_DEPENDENCY_ROOT`
to use a different dependency directory. Meshing and smoothing dependencies
are not needed just to view an existing mesh.

## Generate a mesh

The mesher requires `step_mesher_cpp_matlab`, including its compiled backend,
beside this repository or at `STEPMESHER_ROOT`. This dependency is external.
An example with the settings used for the default rimfree mesh is:

```matlab
addpath('src/dish');
setup_dish;
result = run_dish_mesh('dish_fd023_m', 'order', 4, ...
    'h_rim', 0.03, 'target', [], 'mesh_fraction', 0.10, ...
    'dist_min', 0.02, 'dist_max', 1.0, 'h_max', inf, 'tag', '_rimfree');
[S,fig] = preview_dish(fullfile('src','dish','mesh','out', ...
    'dish_fd023_m_rimfree','dish_fd023_m_rimfree_order4.go3'));
```

All length options use the CAD file's units. The driver defaults use metre
scales; set appropriate sizes explicitly for `_t1` CAD. `h_max` is a reported
diagnostic bound, not a guarantee from the mesher. The driver preserves the
source settings and diagnostics; high-resolution meshes can be large.

## Regenerate CAD

The committed STEP files can be used directly. To rebuild the metre variant,
install the Python requirements and use Gmsh with OpenCASCADE support:

```sh
python3 -m pip install -r src/dish/cad/requirements.txt
python3 src/dish/cad/make_dish_geo.py --f-ratio 0.23 \
  --diameter 6 --thickness 0.01 --units m --out src/dish/cad/dish_fd023_m.geo
gmsh src/dish/cad/dish_fd023_m.geo -0 -o src/dish/cad/dish_fd023_m.step
```

`--units` labels the output; it does not rescale input lengths.
