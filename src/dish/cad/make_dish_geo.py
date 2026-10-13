#!/usr/bin/env python3
"""Emit a gmsh/OpenCASCADE .geo for the HIRAX parabolic reflector.

The reflector in hirax_gum_fid_v0.cst is a zero-thickness swept sheet:

    circle   : radius ddish/2
    profile  : LawX "t", LawZ "t*t/(4*f_length) - f_length", t in [0, ddish/2]
    SweepCurve "reflector", material Aluminum

so the concave reflecting surface is  z = r^2/(4 f) - f,  focus at the origin,
vertex at z = -f.  Here that sheet is given a uniform wall thickness t_wall,
measured along the surface normal, to make a closed solid.

The solid is a revolution of this profile in the (r, z) half plane:

  inner (reflecting)  exact quadratic Bezier.  A parabola IS a quadratic
                      Bezier: control points (0, -f), (R/2, -f), (R, z_rim)
                      reproduce z = r^2/(4f) - f identically.
  rim                 straight segment of length t_wall along the outward
                      surface normal at r = R; revolves to a cone frustum.
  outer               single Bezier of degree 2p+1 for the normal offset.
  axis                segment on r = 0 closing the profile; degenerate under
                      the revolve, so it yields no face.

Why a single Bezier for the offset and not a spline.  The normal offset of a
parabola is not itself a conic, so it has to be approximated.  An interpolating
spline through sampled offset points is accurate but has a knot line per span,
and each knot line is only C2: step_mesher's continuity-aware feature
correction then fragments the outer face along every one of them, which buries
the 2D mesher ("Impossible to recover the edge").  A single Bezier has no
internal knots, so nothing is there to fragment.

Writing the offset in terms of u = r/R makes it polynomial to high accuracy.
With w = u^2 and a = R^2/(4 f^2),

    h(w) = (1 + a w)^(-1/2)
    Z(u) = R^2 u^2/(4 f) - f - t_wall * h(u^2)          even,  degree 2p
    X(u) = R u (1 + (t_wall/(2 f)) * h(u^2))            odd,   degree 2p+1

so only h needs approximating.  It is interpolated at Chebyshev-Lobatto nodes
in w, which include w = 0 and w = 1 and therefore reproduce both endpoints
exactly: the profile meets the axis exactly at z = -f - t_wall, and meets the
rim segment exactly at the normal offset of the rim.  Z even in u also makes
the profile tangent horizontal on the axis, so the revolved pole is smooth.

At the default p = 8 (Bezier degree 17, within OpenCASCADE's degree-25 limit)
the offset profile is within about 1e-6 cm of the exact constant-thickness
offset, i.e. the wall is 1 cm uniform to one part in a million.

The construction is unit-agnostic: use one consistent unit for all lengths.
Defaults are in centimetres, matching the CST model; --units only labels the
output and does not convert the input lengths.
"""

import argparse
import json
import os
from math import comb

import numpy as np
from scipy.integrate import quad


# ---------------------------------------------------------------- exact offset

def offset_point(r, f, t):
    """Normal offset of z = r^2/(4f) by t, away from the concave side."""
    m = r / (2.0 * f)
    s = np.sqrt(1.0 + m * m)
    return r + t * m / s, r * r / (4.0 * f) - f - t / s


# -------------------------------------------------------- Bezier for the offset

def _mono_to_bernstein(c, n):
    """Control points of the degree-n Bezier equal to the polynomial c (ascending)."""
    c = np.concatenate([np.asarray(c, float), np.zeros(n + 1 - len(c))])
    return np.array([sum(comb(i, k) / comb(n, k) * c[k] for k in range(i + 1))
                     for i in range(n + 1)])


def offset_bezier(f, R, t, p):
    """Control points (x, z) of the single Bezier approximating the normal offset."""
    a = R * R / (4.0 * f * f)
    w = (1.0 - np.cos(np.pi * np.arange(p + 1) / p)) / 2.0   # Chebyshev-Lobatto
    cH = np.polyfit(w, (1.0 + a * w) ** -0.5, p)[::-1]       # interpolant in w

    n = 2 * p + 1
    cz = np.zeros(n + 1)
    cz[0] = -f
    cz[2] += R * R / (4.0 * f)
    cx = np.zeros(n + 1)
    cx[1] += R
    for k in range(p + 1):
        cz[2 * k] -= t * cH[k]
        cx[2 * k + 1] += R * (t / (2.0 * f)) * cH[k]

    return _mono_to_bernstein(cx, n), _mono_to_bernstein(cz, n), n


def bezier_eval(b, u):
    """de Casteljau; stable, and well behaved at u = 0 and u = 1."""
    u = np.asarray(u, float)
    c = np.repeat(np.asarray(b, float)[:, None], u.size, axis=1)
    for k in range(len(b) - 1, 0, -1):
        c = (1.0 - u) * c[:k] + u * c[1:k + 1]
    return c[0]


def bezier_deriv(b, u):
    n = len(b) - 1
    d = n * np.diff(b)
    return bezier_eval(d, u)


def offset_error(f, R, t, p, m=200001):
    """Max distance from the Bezier to the exact constant-offset curve.

    Both carry the same parameter u = r/R, so the pointwise distance is an
    upper bound on the geometric deviation.
    """
    bx, bz, _ = offset_bezier(f, R, t, p)
    u = np.linspace(0.0, 1.0, m)
    ex, ez = offset_point(R * u, f, t)
    return float(np.max(np.hypot(bezier_eval(bx, u) - ex, bezier_eval(bz, u) - ez)))


def wall_thickness_range(f, R, t, p, m=2001):
    """Actual wall thickness of the as-built solid, as a sanity check."""
    bx, bz, _ = offset_bezier(f, R, t, p)
    ub = np.linspace(0.0, 1.0, 20001)
    ox, oz = bezier_eval(bx, ub), bezier_eval(bz, ub)
    r = np.linspace(0.0, R, m)
    z = r * r / (4.0 * f) - f
    d = np.array([np.min(np.hypot(ox - ri, oz - zi)) for ri, zi in zip(r, z)])
    return float(d.min()), float(d.max())


# ------------------------------------------------------------- as-built areas

def areas(f, R, t, p):
    """Exact areas of the three faces of the as-built solid."""
    # inner paraboloid, closed form
    a_in = (8.0 * np.pi * f ** 2 / 3.0) * ((1.0 + (R / (2.0 * f)) ** 2) ** 1.5 - 1.0)

    # outer: surface of revolution of the Bezier, 2*pi*x*|C'| du
    bx, bz, _ = offset_bezier(f, R, t, p)

    def dA(u):
        return 2.0 * np.pi * bezier_eval(bx, np.array([u]))[0] * np.hypot(
            bezier_deriv(bx, np.array([u]))[0], bezier_deriv(bz, np.array([u]))[0])

    a_out = quad(dA, 0.0, 1.0, limit=400, epsabs=1e-10, epsrel=1e-13)[0]

    # rim: cone frustum between the inner rim and the Bezier end point
    r1 = R
    z1 = R * R / (4.0 * f) - f
    r2 = float(bezier_eval(bx, np.array([1.0]))[0])
    z2 = float(bezier_eval(bz, np.array([1.0]))[0])
    a_rim = np.pi * (r1 + r2) * np.hypot(r2 - r1, z2 - z1)

    return a_in, a_out, a_rim


# ------------------------------------------------------------------- geo output

def write_geo(path, f, R, t, p, units='cm'):
    bx, bz, n = offset_bezier(f, R, t, p)

    L = []
    L.append('// HIRAX parabolic reflector, solid of revolution.')
    L.append('// Concave reflecting surface: z = r^2/(4*f) - f, focus at the origin,')
    L.append('// exactly the "reflector" sheet of hirax_gum_fid_v0.cst.')
    L.append(f'// f = {f:.10g}, D = {2 * R:.10g}, f/D = {f / (2 * R):.10g} ({units})')
    L.append(f'// wall thickness = {t:.10g}, uniform along the surface normal.')
    L.append(f'// Lengths are in {units}.')
    L.append('SetFactory("OpenCASCADE");')
    L.append('')
    L.append(f'f = {f:.17g};')
    L.append(f'R = {R:.17g};')
    L.append(f't = {t:.17g};')
    L.append('')
    L.append('// Inner (reflecting) profile. A parabola is exactly a quadratic Bezier')
    L.append('// whose middle control point is where the end tangents meet.')
    L.append('Point(1) = {0, 0, -f};')
    L.append('Point(2) = {R/2, 0, -f};')
    L.append('Point(3) = {R, 0, R*R/(4*f) - f};')
    L.append('Bezier(1) = {1, 2, 3};')
    L.append('')
    L.append('// Rim: segment of length t along the outward surface normal at r = R.')
    L.append(f'Point(4) = {{{bx[-1]:.17g}, 0, {bz[-1]:.17g}}};')
    L.append('Line(2) = {3, 4};')
    L.append('')
    L.append(f'// Outer profile: single degree-{n} Bezier for the normal offset, rim to axis.')
    L.append('// One span, no internal knots, so nothing for the mesher to fragment on.')
    for i in range(1, n):
        j = n - i      # written rim-to-axis, i.e. control points reversed
        L.append(f'Point({100 + i}) = {{{bx[j]:.17g}, 0, {bz[j]:.17g}}};')
    L.append('Point(5) = {0, 0, -f-t};')
    ids = ['4'] + [str(100 + i) for i in range(1, n)] + ['5']
    L.append('Bezier(3) = {' + ', '.join(ids) + '};')
    L.append('')
    L.append('// Close the profile along the axis. This edge is degenerate under the')
    L.append('// revolve and contributes no face to the solid.')
    L.append('Line(4) = {5, 1};')
    L.append('')
    L.append('Curve Loop(1) = {1, 2, 3, 4};')
    L.append('Plane Surface(1) = {1};')
    L.append('Extrude { {0, 0, 1}, {0, 0, 0}, 2*Pi } { Surface{1}; }')
    L.append('')
    L.append('// The OCC revolve does not consume its generatrix: Surface(1) survives as')
    L.append('// a loose face, and gmsh would write it into the STEP beside the solid,')
    L.append('// where it inflates the reported CAD area and offers the mesher a')
    L.append('// spurious face to match against. Drop it.')
    L.append('Recursive Delete { Surface{1}; }')
    L.append('')
    with open(path, 'w') as fh:
        fh.write('\n'.join(L) + '\n')

    return n, float(bx[-1]), float(bz[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--f-ratio', type=float, required=True)
    ap.add_argument('--diameter', type=float, default=600.0)
    ap.add_argument('--thickness', type=float, default=1.0)
    ap.add_argument('--units', default='cm',
                    help='label only; the construction is unit-agnostic, so pass '
                         'lengths in whatever unit you want the CAD to carry')
    ap.add_argument('--offset-degree-half', type=int, default=8,
                    help='p; the offset Bezier has degree 2p+1 (default 8 -> 17)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    f = a.f_ratio * a.diameter
    R = a.diameter / 2.0
    t = a.thickness
    p = a.offset_degree_half

    n, r_rim_out, z_rim_out = write_geo(a.out, f, R, t, p, a.units)
    a_in, a_out, a_rim = areas(f, R, t, p)
    tmin, tmax = wall_thickness_range(f, R, t, p)

    print(f'wrote {a.out}')
    print(f'  f/D             {a.f_ratio:.10g}')
    print(f'  f               {f:.10g} {a.units}')
    print(f'  D               {a.diameter:.10g} {a.units}')
    print(f'  wall thickness  {t:.10g} {a.units} along the normal')
    print(f'  vertex z        {-f:.10g} (inner), {-f - t:.10g} (outer)')
    print(f'  rim inner       r = {R:.10g}, z = {R * R / (4 * f) - f:.10g}')
    print(f'  rim outer       r = {r_rim_out:.10g}, z = {z_rim_out:.10g}')
    print(f'  offset Bezier   degree {n}, max deviation from the exact offset '
          f'{offset_error(f, R, t, p):.3e} {a.units}')
    print(f'  as-built wall   {tmin:.9f} to {tmax:.9f} {a.units}')
    print(f'  exact areas     inner {a_in:.6f}')
    print(f'                  outer {a_out:.6f}')
    print(f'                  rim   {a_rim:.6f}')
    print(f'                  TOTAL {a_in + a_out + a_rim:.6f} {a.units}^2')

    # Companion record, so the meshing driver can validate against exact values
    # instead of a CAD-reported area.
    meta = {
        'step_stem': os.path.splitext(os.path.basename(a.out))[0],
        'units': a.units,
        'f_ratio': a.f_ratio,
        'focal_length': f,
        'diameter': a.diameter,
        'aperture_radius': R,
        'wall_thickness': t,
        'wall_thickness_min': tmin,
        'wall_thickness_max': tmax,
        'offset_bezier_degree': n,
        'offset_bezier_max_deviation': offset_error(f, R, t, p),
        'vertex_z_inner': -f,
        'vertex_z_outer': -f - t,
        'rim_inner_r': R,
        'rim_inner_z': R * R / (4 * f) - f,
        'rim_outer_r': r_rim_out,
        'rim_outer_z': z_rim_out,
        'area_inner_exact': a_in,
        'area_outer_exact': a_out,
        'area_rim_exact': a_rim,
        'area_total_exact': a_in + a_out + a_rim,
    }
    jpath = os.path.splitext(a.out)[0] + '.json'
    with open(jpath, 'w') as fh:
        json.dump(meta, fh, indent=2, sort_keys=True)
    print(f'  wrote {jpath}')


if __name__ == '__main__':
    main()
