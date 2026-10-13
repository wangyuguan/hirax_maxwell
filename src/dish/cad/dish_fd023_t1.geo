// HIRAX parabolic reflector, solid of revolution.
// Concave reflecting surface: z = r^2/(4*f) - f, focus at the origin,
// exactly the "reflector" sheet of hirax_gum_fid_v0.cst.
// f = 138 cm, D = 600 cm, f/D = 0.23
// wall thickness = 1 cm, uniform along the surface normal.
// Lengths in centimetres, as in the CST model.
SetFactory("OpenCASCADE");

f = 138;
R = 300;
t = 1;

// Inner (reflecting) profile. A parabola is exactly a quadratic Bezier
// whose middle control point is where the end tangents meet.
Point(1) = {0, 0, -f};
Point(2) = {R/2, 0, -f};
Point(3) = {R, 0, R*R/(4*f) - f};
Bezier(1) = {1, 2, 3};

// Rim: segment of length t along the outward surface normal at r = R.
Point(4) = {300.73593101176175, 0, 24.366421730048682};
Line(2) = {3, 4};

// Outer profile: single degree-17 Bezier for the normal offset, rim to axis.
// One span, no internal knots, so nothing for the mesher to fragment on.
Point(101) = {283.0690271353331, 0, 5.1632666462125139};
Point(102) = {265.40011509993161, 0, -12.841888002270176};
Point(103) = {247.72892865536497, 0, -29.648966999123751};
Point(104) = {230.0553123388577, 0, -45.257975562131186};
Point(105) = {212.3788589954745, 0, -59.668696798100832};
Point(106) = {194.69939667182717, 0, -72.881019308735844};
Point(107) = {177.01661983967313, 0, -84.894727312608012};
Point(108) = {159.33023230637801, 0, -95.709478148931098};
Point(109) = {141.64001805816659, 0, -105.32485539931965};
Point(110) = {123.94583385788425, 0, -113.74035671409911};
Point(111) = {106.24765072594272, 0, -120.95536475808994};
Point(112) = {88.545636517101045, 0, -126.96917188287418};
Point(113) = {70.840212745594741, 0, -131.78106364569683};
Point(114) = {53.132048071309178, 0, -135.39042195184152};
Point(115) = {35.421994884910482, 0, -137.7968073172805};
Point(116) = {17.710997442455241, 0, -139};
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

