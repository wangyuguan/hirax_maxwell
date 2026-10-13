// HIRAX parabolic reflector, solid of revolution.
// Concave reflecting surface: z = r^2/(4*f) - f, focus at the origin,
// exactly the "reflector" sheet of hirax_gum_fid_v0.cst.
// f = 1.26, D = 6, f/D = 0.21 (m)
// wall thickness = 0.01, uniform along the surface normal.
// Lengths are in m.
SetFactory("OpenCASCADE");

f = 1.26;
R = 3;
t = 0.01;

// Inner (reflecting) profile. A parabola is exactly a quadratic Bezier
// whose middle control point is where the end tangents meet.
Point(1) = {0, 0, -f};
Point(2) = {R/2, 0, -f};
Point(3) = {R, 0, R*R/(4*f) - f};
Bezier(1) = {1, 2, 3};

// Rim: segment of length t along the outward surface normal at r = R.
Point(4) = {3.0076570486478955, 0, 0.51928236485005308};
Line(2) = {3, 4};

// Outer profile: single degree-17 Bezier for the normal offset, rim to axis.
// One span, no internal knots, so nothing for the mesher to fragment on.
Point(101) = {2.8310001060042005, 0, 0.30897652065353387};
Point(102) = {2.6543228913756751, 0, 0.11179023070222691};
Point(103) = {2.4776216575910546, 0, -0.072275496951128823};
Point(104) = {2.300896288447694, 0, -0.24322241625237037};
Point(105) = {2.1241393818922236, 0, -0.40104693532479951};
Point(106) = {1.9473505261232327, 0, -0.54574911679203075};
Point(107) = {1.7705253983112876, 0, -0.67732736218953982};
Point(108) = {1.5936593182867431, 0, -0.79577781939579717};
Point(109) = {1.4167489881863873, 0, -0.90109591861390004};
Point(110) = {1.2397917363820021, 0, -0.99327620203894162};
Point(111) = {1.0627857840531199, 0, -1.0723113691162587};
Point(112) = {0.88573173614339873, 0, -1.138192254561293};
Point(113) = {0.7086338526364695, 0, -1.1909090353467753};
Point(114) = {0.53150019985379526, 0, -1.2304529383258875};
Point(115) = {0.35434173669467783, 0, -1.2568176461086291};
Point(116) = {0.17717086834733892, 0, -1.27};
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

