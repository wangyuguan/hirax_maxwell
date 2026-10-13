// HIRAX parabolic reflector, solid of revolution.
// Concave reflecting surface: z = r^2/(4*f) - f, focus at the origin,
// exactly the "reflector" sheet of hirax_gum_fid_v0.cst.
// f = 1.38, D = 6, f/D = 0.23 (m)
// wall thickness = 0.01, uniform along the surface normal.
// Lengths are in m.
SetFactory("OpenCASCADE");

f = 1.3800000000000001;
R = 3;
t = 0.01;

// Inner (reflecting) profile. A parabola is exactly a quadratic Bezier
// whose middle control point is where the end tangents meet.
Point(1) = {0, 0, -f};
Point(2) = {R/2, 0, -f};
Point(3) = {R, 0, R*R/(4*f) - f};
Bezier(1) = {1, 2, 3};

// Rim: segment of length t along the outward surface normal at r = R.
Point(4) = {3.0073593101176179, 0, 0.24366421730048654};
Line(2) = {3, 4};

// Outer profile: single degree-17 Bezier for the normal offset, rim to axis.
// One span, no internal knots, so nothing for the mesher to fragment on.
Point(101) = {2.8306902713533311, 0, 0.051632666462124957};
Point(102) = {2.6540011509993153, 0, -0.12841888002270185};
Point(103) = {2.47728928655365, 0, -0.29648966999123771};
Point(104) = {2.300553123388577, 0, -0.45257975562131209};
Point(105) = {2.1237885899547444, 0, -0.59668696798100851};
Point(106) = {1.9469939667182721, 0, -0.72881019308735873};
Point(107) = {1.7701661983967312, 0, -0.84894727312608043};
Point(108) = {1.5933023230637797, 0, -0.95709478148931126};
Point(109) = {1.4164001805816662, 0, -1.0532485539931966};
Point(110) = {1.2394583385788427, 0, -1.1374035671409912};
Point(111) = {1.0624765072594271, 0, -1.2095536475808997};
Point(112) = {0.88545636517101045, 0, -1.269691718828742};
Point(113) = {0.70840212745594744, 0, -1.3178106364569686};
Point(114) = {0.53132048071309168, 0, -1.3539042195184152};
Point(115) = {0.35421994884910485, 0, -1.3779680731728052};
Point(116) = {0.17710997442455242, 0, -1.3900000000000001};
Point(5) = {0, 0, -f-t};
Bezier(3) = {4, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 5};

// Close the profile along the axis. This edge is degenerate under the
// revolve and contributes no face to the solid.
Line(4) = {5, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Extrude { {0, 0, 1}, {0, 0, 0}, 2*Pi } { Surface{1}; }

// The OCC revolve does not consume its generatrix: Surface(1) survives as
// a loose face, and gmsh would write it into the STEP beside the solid,
// where it inflates the reported CAD area and offers the mesher a
// spurious face to match against. Drop it.
Recursive Delete { Surface{1}; }

