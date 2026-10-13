// HIRAX parabolic reflector, solid of revolution.
// Concave reflecting surface: z = r^2/(4*f) - f, focus at the origin,
// exactly the "reflector" sheet of hirax_gum_fid_v0.cst.
// f = 126 cm, D = 600 cm, f/D = 0.21
// wall thickness = 1 cm, uniform along the surface normal.
// Lengths in centimetres, as in the CST model.
SetFactory("OpenCASCADE");

f = 126;
R = 300;
t = 1;

// Inner (reflecting) profile. A parabola is exactly a quadratic Bezier
// whose middle control point is where the end tangents meet.
Point(1) = {0, 0, -f};
Point(2) = {R/2, 0, -f};
Point(3) = {R, 0, R*R/(4*f) - f};
Bezier(1) = {1, 2, 3};

// Rim: segment of length t along the outward surface normal at r = R.
Point(4) = {300.76570486478965, 0, 51.9282364850053};
Line(2) = {3, 4};

// Outer profile: single degree-17 Bezier for the normal offset, rim to axis.
// One span, no internal knots, so nothing for the mesher to fragment on.
Point(101) = {283.10001060042009, 0, 30.897652065353359};
Point(102) = {265.43228913756758, 0, 11.179023070222691};
Point(103) = {247.76216575910544, 0, -7.2275496951128808};
Point(104) = {230.08962884476944, 0, -24.322241625237034};
Point(105) = {212.41393818922239, 0, -40.104693532479956};
Point(106) = {194.73505261232327, 0, -54.574911679203069};
Point(107) = {177.05253983112877, 0, -67.732736218953974};
Point(108) = {159.36593182867432, 0, -79.577781939579737};
Point(109) = {141.67489881863875, 0, -90.10959186139003};
Point(110) = {123.97917363820022, 0, -99.327620203894156};
Point(111) = {106.278578405312, 0, -107.23113691162588};
Point(112) = {88.573173614339879, 0, -113.8192254561293};
Point(113) = {70.863385263646961, 0, -119.09090353467752};
Point(114) = {53.150019985379529, 0, -123.04529383258875};
Point(115) = {35.434173669467789, 0, -125.68176461086291};
Point(116) = {17.717086834733895, 0, -127};
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

