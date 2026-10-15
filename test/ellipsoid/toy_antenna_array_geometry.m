function [model,physical] = toy_antenna_array_geometry(norder)
%TOY_ANTENNA_ARRAY_GEOMETRY Cheap closed surfaces for nested-array experiments.
%   [model,physical] = toy_antenna_array_geometry(4)
%   model = toy_antenna_array_geometry(0) returns analytic parameters only.
%
%   Lengths are in wavelengths, lambda=1 and k=2*pi. One station contains
%   an oblate ellipsoid and four smaller radial ellipsoids above it. The
%   inner sphere encloses the four leaves but excludes the body. An outer
%   oblate ellipsoid encloses the complete station and inner sphere. Four
%   translations form a 2-by-2 array with disjoint enclosing surfaces.
%   The body-ground gap/diameter ratio is 1/6, using the provisional 1 m
%   clearance for a 6 m dish. This preserves geometric, not electrical,
%   clearance at the real operating frequency; this toy keeps kD=2*pi.
%
%   physical.body, physical.leaves, physical.feed, and physical.station are
%   surfer objects for ONE station. Reuse this station by translations from
%   model.station_shifts. All ellipsoids use six tensor-product patches.
%   norder=4 is an initial mesh, not a certified scattering accuracy.
%   This model tests compression/coupling, not paraboloid focusing or ports.
%   The plane z=0 is a geometric reference; no ground response is included.

    if nargin<1 || isempty(norder), norder = 4; end
    validateattributes(norder,{'numeric'},{'scalar','integer','nonnegative'});
    model.wavelength = 1;
    model.zk = 2*pi;
    model.reference_dish_diameter_m = 6;
    model.reference_ground_clearance_m = 1; % Provisional value from Ben.
    model.ground_z = 0;
    model.body_axes = [0.50;0.50;0.16];
    model.geometric_scale_m = model.reference_dish_diameter_m/(2*model.body_axes(1));
    model.ground_clearance = model.reference_ground_clearance_m/model.geometric_scale_m;
    model.body_center = [0;0;model.ground_z+model.ground_clearance+model.body_axes(3)];
    model.leaf_axes = [0.085;0.045;0.035];
    model.leaf_radial_offset = 0.15;
    model.leaf_angles = pi/4+(0:3)*pi/2;
    model.inner_center = model.body_center+[0;0;0.55];
    model.leaf_centers = [model.leaf_radial_offset*cos(model.leaf_angles); ...
        model.leaf_radial_offset*sin(model.leaf_angles); ...
        model.inner_center(3)*ones(1,4)];
    model.inner_radius = 0.28;
    model.outer_center = model.body_center+[0;0;0.29];
    model.outer_axes = [0.80;0.80;0.60];
    model.spacing = 1.90;
    model.station_shifts = model.spacing* ...
        [-0.5,0.5,-0.5,0.5;-0.5,-0.5,0.5,0.5;0,0,0,0];
    model.ground_response_included = false;
    model.norder = norder;
    model.description = ['Four small ellipsoids coupled to one oblate ellipsoid; ', ...
        'inner sphere and outer ellipsoid; four translated stations; h/D=1/6; lambda=1'];

    % Conservative analytic separation checks apply to the entire surfaces,
    % independently of the plotted samples or the chosen discretization.
    model.clearances.body_to_ground = model.body_center(3)-model.body_axes(3)-model.ground_z;
    model.clearances.inner_to_body = model.inner_center(3)-model.inner_radius ...
        -model.body_center(3)-model.body_axes(3);
    model.clearances.outer_to_ground = model.outer_center(3)-model.outer_axes(3)-model.ground_z;
    model.clearances.outer_to_neighbor = model.spacing-2*max(model.outer_axes(1:2));
    % The coaxial body's bounding cylinder is inside the outer ellipsoid.
    % The sphere uses a triangle inequality in scaled outer coordinates.
    % Multiplication by the shortest semiaxis gives distance lower bounds.
    body_bound = hypot(max(model.body_axes(1:2))/min(model.outer_axes(1:2)), ...
        (abs(model.body_center(3)-model.outer_center(3))+model.body_axes(3))/model.outer_axes(3));
    inner_bound = norm((model.inner_center-model.outer_center)./model.outer_axes) ...
        +model.inner_radius/min(model.outer_axes);
    model.clearances.body_inside_outer = min(model.outer_axes)*(1-body_bound);
    model.clearances.inner_inside_outer = min(model.outer_axes)*(1-inner_bound);
    model.clearances.leaves_inside_inner = model.inner_radius ...
        -model.leaf_radial_offset-max(model.leaf_axes);
    model.clearances.between_leaves = sqrt(2)*model.leaf_radial_offset ...
        -2*max(model.leaf_axes);
    assert(all(cell2mat(struct2cell(model.clearances))>0));
    physical = [];
    if norder==0, return; end

    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    run(fullfile(project_root,'startup.m'));
    physical = struct();
    physical.body = geometries.ellipsoid(model.body_axes,[1;1;1], ...
        model.body_center,norder,11);
    base = geometries.ellipsoid(model.leaf_axes,[1;1;1],[0;0;0],norder,11);
    physical.leaves = cell(1,4);
    for j = 1:4
      a = model.leaf_angles(j);
      R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
      physical.leaves{j} = affine_transf(base,R,model.leaf_centers(:,j));
    end
    physical.feed = merge([physical.leaves{:}]);
    physical.station = merge([physical.body,physical.feed]);
    model.stats.patches_per_component = base.npatches;
    model.stats.nodes_per_component = base.npts;
    model.stats.station_nodes = physical.station.npts;
    model.stats.station_nrccie_dofs = 3*physical.station.npts;
    model.stats.array_nodes = 4*physical.station.npts;
    model.stats.array_nrccie_dofs = 12*physical.station.npts;
    assert(base.npatches==6 && base.npts==6*(norder+1)^2);
    assert(physical.station.npts==5*base.npts);
    assert(all(isfinite(physical.station.r(:))));
    assert(all(physical.station.wts>0));
    assert(max(abs(vecnorm(physical.station.n,2,1)-1))<1e-10);
    fprintf('Toy geometry p%d: %d nodes/component, %d nodes/station, %d array nodes\n', ...
        norder,base.npts,physical.station.npts,model.stats.array_nodes);
    fprintf('Uncompressed NRCCIE DOFs: %d per station; %d for four stations\n', ...
        model.stats.station_nrccie_dofs,model.stats.array_nrccie_dofs);
end
