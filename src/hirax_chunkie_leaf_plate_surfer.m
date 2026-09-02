function [S,parts] = hirax_chunkie_leaf_plate_surfer(thickness,opts)
%HIRAX_CHUNKIE_LEAF_PLATE_SURFER Smooth HIRAX-shaped thin plate.
%
% The physical geometry is generated directly in the original HIRAX length
% scale.  The planform is translated to the origin, but it is not divided by
% the nominal leaf radius.  No saved msh, rquad, or surfer mesh is loaded.
% Chunkie's smoother rounds every vertex of the input polygon, including the
% two lower corners, and produces one fixed master curve.  Changing norder
% only resamples that same master curve.

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'norder') || isempty(opts.norder)
    opts.norder = 4;
end
if ~isfield(opts,'chunkie_order') || isempty(opts.chunkie_order)
    opts.chunkie_order = 20;
end
if ~isfield(opts,'chunkie_n0') || isempty(opts.chunkie_n0)
    opts.chunkie_n0 = 3;
end
if ~isfield(opts,'chunkie_nchs') || isempty(opts.chunkie_nchs)
    opts.chunkie_nchs = 3;
end
if ~isfield(opts,'chunkie_newton_iterations') || ...
        isempty(opts.chunkie_newton_iterations)
    opts.chunkie_newton_iterations = 30;
end
if ~isfield(opts,'rim_width') || isempty(opts.rim_width)
    opts.rim_width = 0.028128271246*69.25;
end
if ~isfield(opts,'cap_collar_width') || isempty(opts.cap_collar_width)
    opts.cap_collar_width = opts.rim_width;
end
if ~isfield(opts,'cap_mesh_spacing') || isempty(opts.cap_mesh_spacing)
    opts.cap_mesh_spacing = 0.09*69.25;
end
if ~isfield(opts,'cap_mesh_spacing_center') || ...
        isempty(opts.cap_mesh_spacing_center)
    opts.cap_mesh_spacing_center = opts.cap_mesh_spacing;
end
if ~isfield(opts,'cap_mesh_spacing_side') || ...
        isempty(opts.cap_mesh_spacing_side)
    opts.cap_mesh_spacing_side = opts.cap_mesh_spacing;
end
if ~isfield(opts,'cap_mesh_side_start') || ...
        isempty(opts.cap_mesh_side_start)
    opts.cap_mesh_side_start = 0.55;
end
if ~isfield(opts,'wall_profile_refinement') || ...
        isempty(opts.wall_profile_refinement)
    opts.wall_profile_refinement = 1;
end
if ~isfield(opts,'outline_refinement') || ...
        isempty(opts.outline_refinement)
    opts.outline_refinement = 1;
end
validateattributes(thickness,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.norder,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',4});
validateattributes(opts.chunkie_order,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',10});
validateattributes(opts.chunkie_n0,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',2});
validateattributes(opts.chunkie_nchs,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1});
validateattributes(opts.chunkie_newton_iterations,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1});
validateattributes(opts.rim_width,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.cap_collar_width,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.cap_mesh_spacing,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.cap_mesh_spacing_center,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.cap_mesh_spacing_side,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(opts.cap_mesh_side_start,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',1});
validateattributes(opts.wall_profile_refinement,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1});
validateattributes(opts.outline_refinement,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1});
norder = opts.norder;
[outline,design] = fixed_chunkie_outline(opts);
number_of_master_outline_panels = outline.number_of_panels;
if opts.outline_refinement>1
    outline_order = outline.order;
    [legendre_nodes,~,values_to_coefficients] = ...
        lege.exps(outline_order);
    local_q = 0.5*(legendre_nodes(:).'+1);
    number_of_refined_panels = opts.outline_refinement* ...
        number_of_master_outline_panels;
    refined_position_coefficients = zeros( ...
        2,outline_order,number_of_refined_panels);

    refined_panel = 0;
    for master_panel = 1:number_of_master_outline_panels
        for child = 1:opts.outline_refinement
            refined_panel = refined_panel+1;
            parent_q = (child-1+local_q)/ ...
                opts.outline_refinement;
            position = chunkie_panel_values( ...
                outline,master_panel,parent_q);
            refined_position_coefficients(:,:,refined_panel) = ...
                (values_to_coefficients*position.').';
        end
    end

    outline.number_of_panels = number_of_refined_panels;
    outline.position_coefficients = refined_position_coefficients;
    [outline.panel_starts,outline.panel_start_derivatives] = ...
        chunkie_panel_starts(outline);
end
number_of_outline_panels = outline.number_of_panels;
outline_points = outline.panel_starts;

outline_center = [0;0];
mean_outline_radius = mean(vecnorm( ...
    outline_points-outline_center,2,1));
inner_scale = 1-opts.rim_width/mean_outline_radius;
core_scale = inner_scale-opts.cap_collar_width/mean_outline_radius;

inner_points = outline_center+inner_scale*( ...
    outline_points-outline_center);
core_vertices = outline_center+core_scale*( ...
    outline_points-outline_center);

[core_nodes,core_faces] = triangulate_core(core_vertices, ...
    opts.cap_mesh_spacing_center,opts.cap_mesh_spacing_side, ...
    opts.cap_mesh_side_start);
number_of_core_triangles = size(core_faces,2);
number_of_cap_collar_quads = number_of_outline_panels;
number_of_wall_profile_panels = opts.wall_profile_refinement;
number_of_wall_quads_per_half = number_of_outline_panels* ...
    number_of_wall_profile_panels;

triangle_uv = koorn.rv_nodes(norder);
triangle_u = triangle_uv(1,:);
triangle_v = triangle_uv(2,:);
number_of_triangle_nodes = size(triangle_uv,2);

quad_uv = polytens.lege.nodes(norder);
quad_u = 0.5*(quad_uv(1,:)+1);
quad_v = 0.5*(quad_uv(2,:)+1);
number_of_quad_nodes = size(quad_uv,2);
half_thickness = thickness/2;

top_core_values = zeros(12,number_of_core_triangles* ...
    number_of_triangle_nodes);
bottom_core_values = zeros(size(top_core_values));
for triangle = 1:number_of_core_triangles
    ids = (triangle-1)*number_of_triangle_nodes+ ...
        (1:number_of_triangle_nodes);
    vertices = core_nodes(:,core_faces(:,triangle));
    top_core_values(:,ids) = planar_triangle_values( ...
        vertices,half_thickness,triangle_u,triangle_v);
    bottom_core_values(:,ids) = planar_triangle_values( ...
        vertices(:,[1 3 2]),-half_thickness,triangle_u,triangle_v);
end

top_collar_values = zeros(12,number_of_cap_collar_quads* ...
    number_of_quad_nodes);
bottom_collar_values = zeros(size(top_collar_values));
upper_wall_values = zeros(12,number_of_wall_quads_per_half* ...
    number_of_quad_nodes);
lower_wall_values = zeros(size(upper_wall_values));

for panel = 1:number_of_outline_panels
    next_panel = mod(panel,number_of_outline_panels)+1;
    collar_ids = (panel-1)*number_of_quad_nodes+ ...
        (1:number_of_quad_nodes);

    [top_position,top_du,top_dv] = cap_collar_values( ...
        outline,inner_scale,core_vertices,panel,next_panel, ...
        quad_u,quad_v,half_thickness,+1);
    top_normal = normalized_cross(top_du,top_dv);
    top_collar_values(:,collar_ids) = ...
        [top_position;top_du;top_dv;top_normal];

    [bottom_position,bottom_du,bottom_dv] = cap_collar_values( ...
        outline,inner_scale,core_vertices,panel,next_panel, ...
        quad_u,quad_v,-half_thickness,-1);
    bottom_normal = normalized_cross(bottom_du,bottom_dv);
    bottom_collar_values(:,collar_ids) = ...
        [bottom_position;bottom_du;bottom_dv;bottom_normal];

    for profile_panel = 1:number_of_wall_profile_panels
        wall_patch = (panel-1)*number_of_wall_profile_panels+ ...
            profile_panel;
        wall_ids = (wall_patch-1)*number_of_quad_nodes+ ...
            (1:number_of_quad_nodes);
        profile_interval = (profile_panel-1:profile_panel)/ ...
            number_of_wall_profile_panels;

        [upper_position,upper_du,upper_dv] = wall_values( ...
            outline,inner_scale,panel,quad_u,quad_v, ...
            half_thickness,+1,profile_interval);
        upper_normal = normalized_cross(upper_du,upper_dv);
        upper_wall_values(:,wall_ids) = [upper_position;upper_du; ...
            upper_dv;upper_normal];

        [lower_position,lower_du,lower_dv] = wall_values( ...
            outline,inner_scale,panel,quad_u,quad_v, ...
            half_thickness,-1,profile_interval);
        lower_normal = normalized_cross(lower_du,lower_dv);
        lower_wall_values(:,wall_ids) = [lower_position;lower_du; ...
            lower_dv;lower_normal];
    end
end

top_core = surfer(number_of_core_triangles,norder,top_core_values,1);
top_collar = surfer(number_of_cap_collar_quads,norder, ...
    top_collar_values,11);
upper_wall = surfer(number_of_wall_quads_per_half,norder, ...
    upper_wall_values,11);
lower_wall = surfer(number_of_wall_quads_per_half,norder, ...
    lower_wall_values,11);
bottom_collar = surfer(number_of_cap_collar_quads,norder, ...
    bottom_collar_values,11);
bottom_core = surfer(number_of_core_triangles,norder, ...
    bottom_core_values,1);

S = merge([top_core,top_collar,upper_wall,lower_wall, ...
    bottom_collar,bottom_core]);

parts = struct();
offset = 0;
parts.top_core = offset+(1:top_core.npatches);
offset = parts.top_core(end);
parts.top_collar = offset+(1:top_collar.npatches);
offset = parts.top_collar(end);
parts.upper_wall = offset+(1:upper_wall.npatches);
offset = parts.upper_wall(end);
parts.lower_wall = offset+(1:lower_wall.npatches);
offset = parts.lower_wall(end);
parts.bottom_collar = offset+(1:bottom_collar.npatches);
offset = parts.bottom_collar(end);
parts.bottom_core = offset+(1:bottom_core.npatches);
parts.top = [parts.top_core parts.top_collar];
parts.wall = [parts.upper_wall parts.lower_wall];
parts.upper_wall_by_profile = reshape(parts.upper_wall, ...
    number_of_wall_profile_panels,number_of_outline_panels);
parts.lower_wall_by_profile = reshape(parts.lower_wall, ...
    number_of_wall_profile_panels,number_of_outline_panels);
parts.bottom = [parts.bottom_collar parts.bottom_core];
parts.outline_points = outline_points;
parts.outline_derivatives = outline.panel_start_derivatives;
parts.outline_position_coefficients = outline.position_coefficients;
parts.chunkie_order = outline.order;
parts.chunkie_raw_vertices = outline.raw_vertices;
parts.inner_points = inner_points;
parts.core_vertices = core_vertices;
parts.core_nodes = core_nodes;
parts.core_faces = core_faces;
parts.number_of_outline_panels = number_of_outline_panels;
parts.number_of_master_outline_panels = ...
    number_of_master_outline_panels;
parts.outline_refinement = opts.outline_refinement;
parts.number_of_wall_profile_panels = number_of_wall_profile_panels;
parts.number_of_core_triangles = number_of_core_triangles;
parts.number_of_cap_collar_layers = 1;
parts.rim_width = opts.rim_width;
parts.cap_collar_width = opts.cap_collar_width;
parts.thickness = thickness;
parts.inner_scale = inner_scale;
parts.core_scale = core_scale;
parts.design = design;
parts.leaf_radius = design.leaf_scale;
parts.profile = 'fixed_quintic_c2';
parts.outline = 'fixed_chunkie_smoothed_no_corners';
parts.options = opts;

parts.max_position_mismatch = shared_edge_position_error(S,parts,norder);
parts.max_normal_mismatch_degrees = shared_edge_normal_error(S,parts,norder);

end


function [outline,design] = fixed_chunkie_outline(opts)
% Cache one deterministic master curve inside this MATLAB process.  The
% cache avoids rerunning Chunkie's Newton smoother at every surface order.
persistent cached_parameters cached_outline cached_design
parameters = [opts.chunkie_order opts.chunkie_n0 opts.chunkie_nchs ...
    opts.chunkie_newton_iterations];
if ~isempty(cached_parameters) && isequal(parameters,cached_parameters)
    outline = cached_outline;
    design = cached_design;
    return
end

% Original single-leaf dimensions, in the same units as the HIRAX script.
x0 = 0;
y0 = 80.35;
x1 = -49.25;
y1 = 80.35;
x2 = 49.25;
y2 = 80.35;
x3 = 0;
y3 = 2.83;
circle_radius = 20.0;
bottom_half_width = 13.8;
ellipse_radius_x = 69.25;
ellipse_radius_y = 54.41;

bottom_left = [x3-bottom_half_width;y3+bottom_half_width];
bottom_right = [x3+bottom_half_width;y3+bottom_half_width];
right_circle_start = [x2+circle_radius*cos(-pi/4); ...
    y2+circle_radius*sin(-pi/4)];
right_ellipse_join = [x2+circle_radius;y2];
left_ellipse_join = [x1-circle_radius;y1];
left_circle_end = [x1+circle_radius*cos(5*pi/4); ...
    y1+circle_radius*sin(5*pi/4)];

n0 = opts.chunkie_n0;
bottom = open_line(bottom_left,bottom_right,n0);
right_line = open_line(bottom_right,right_circle_start,ceil(7*n0/4));
right_circle = open_circle([x2;y2],circle_radius,-pi/4,0,n0);
upper_ellipse = open_ellipse([x0;y0],ellipse_radius_x, ...
    ellipse_radius_y,0,pi,ceil(30*n0/4));
left_circle = open_circle([x1;y1],circle_radius,pi,5*pi/4,n0);
left_line = open_line(left_circle_end,bottom_left,ceil(7*n0/4));
raw_vertices = [bottom right_line right_circle upper_ellipse ...
    left_circle left_line];

smoother_options = struct();
smoother_options.k = opts.chunkie_order;
smoother_options.n_newton = opts.chunkie_newton_iterations;
smoother_options.nchs = opts.chunkie_nchs;
master = chnk.smoother.smooth(raw_vertices,smoother_options);
master = sort(master);

leaf_center = [0;0.5*(bottom_left(2)+y0+ellipse_radius_y)];
leaf_scale = ellipse_radius_x;
number_of_panels = master.nch;
[~,~,values_to_coefficients] = lege.exps(master.k);
position_coefficients = zeros(2,master.k,number_of_panels);
for panel = 1:number_of_panels
    values = squeeze(master.r(:,:,panel));
    values = values-leaf_center;
    position_coefficients(:,:,panel) = ...
        (values_to_coefficients*values.').';
end

outline = struct();
outline.order = master.k;
outline.number_of_panels = number_of_panels;
outline.position_coefficients = position_coefficients;
outline.raw_vertices = raw_vertices-leaf_center;
[outline.panel_starts,outline.panel_start_derivatives] = ...
    chunkie_panel_starts(outline);

design = struct();
design.bottom_left = bottom_left-leaf_center;
design.bottom_right = bottom_right-leaf_center;
design.right_circle_start = right_circle_start-leaf_center;
design.right_ellipse_join = right_ellipse_join-leaf_center;
design.left_ellipse_join = left_ellipse_join-leaf_center;
design.left_circle_end = left_circle_end-leaf_center;
design.leaf_center_raw = leaf_center;
design.leaf_scale = leaf_scale;

cached_parameters = parameters;
cached_outline = outline;
cached_design = design;
end


function points = open_line(point_a,point_b,number_of_points)
q = (0:number_of_points-1)/number_of_points;
points = point_a+(point_b-point_a).*q;
end


function points = open_circle(center,radius,theta_start,theta_end, ...
        number_of_points)
q = (0:number_of_points-1)/number_of_points;
theta = theta_start+(theta_end-theta_start)*q;
points = center+radius*[cos(theta);sin(theta)];
end


function points = open_ellipse(center,radius_x,radius_y,theta_start, ...
        theta_end,number_of_points)
q = (0:number_of_points-1)/number_of_points;
theta = theta_start+(theta_end-theta_start)*q;
points = center+[radius_x*cos(theta);radius_y*sin(theta)];
end


function [starts,derivatives] = chunkie_panel_starts(outline)
starts = zeros(2,outline.number_of_panels);
derivatives = zeros(size(starts));
for panel = 1:outline.number_of_panels
    [starts(:,panel),derivatives(:,panel)] = ...
        chunkie_panel_values(outline,panel,0);
end
end


function [position,derivative_u,derivative_v] = cap_collar_values( ...
        outline,inner_scale,core_vertices, ...
        panel,next_panel,u,v,height,orientation)
if orientation>0
    panel_parameter = u;
    collar_parameter = v;
else
    collar_parameter = u;
    panel_parameter = v;
end

[outer,outer_derivative] = chunkie_panel_values( ...
    outline,panel,panel_parameter);
inner = inner_scale*outer;
inner_derivative = inner_scale*outer_derivative;
core_start = core_vertices(:,panel);
core_end = core_vertices(:,next_panel);
core = core_start+(core_end-core_start).*panel_parameter;
core_derivative = repmat(core_end-core_start,1,numel(panel_parameter));

xy = (1-collar_parameter).*inner+collar_parameter.*core;
dxy_dpanel = (1-collar_parameter).*inner_derivative+ ...
    collar_parameter.*core_derivative;
dxy_dcollar = core-inner;
position = [xy;height*ones(1,numel(panel_parameter))];
if orientation>0
    derivative_u = 0.5*[dxy_dpanel;zeros(1,numel(panel_parameter))];
    derivative_v = 0.5*[dxy_dcollar;zeros(1,numel(panel_parameter))];
else
    derivative_u = 0.5*[dxy_dcollar;zeros(1,numel(panel_parameter))];
    derivative_v = 0.5*[dxy_dpanel;zeros(1,numel(panel_parameter))];
end
end


function [position,derivative_u,derivative_v] = wall_values( ...
        outline,inner_scale,panel,u,v, ...
        half_thickness,half_id,profile_interval)
[outer,outer_derivative] = chunkie_panel_values(outline,panel,v);
inner = inner_scale*outer;
inner_derivative = inner_scale*outer_derivative;

profile_span = profile_interval(2)-profile_interval(1);
parent_u = profile_interval(1)+profile_span*u;
if half_id>0
    profile_parameter = parent_u;
    profile_derivative_scale = profile_span;
    profile_sign = 1;
else
    profile_parameter = 1-parent_u;
    profile_derivative_scale = -profile_span;
    profile_sign = -1;
end
[radial,height,radial_derivative,height_derivative] = ...
    fixed_quintic_profile(profile_parameter);
radial_derivative = profile_derivative_scale*radial_derivative;
height_derivative = profile_derivative_scale*height_derivative;

xy = inner+radial.*(outer-inner);
dxy_ds = radial_derivative.*(outer-inner);
dxy_dt = inner_derivative+radial.*( ...
    outer_derivative-inner_derivative);
z = profile_sign*half_thickness*height;
dz_ds = profile_sign*half_thickness*height_derivative;
position = [xy;z];
derivative_u = 0.5*[dxy_ds;dz_ds];
derivative_v = 0.5*[dxy_dt;zeros(1,numel(v))];
end


function [radial,height,radial_derivative,height_derivative] = ...
        fixed_quintic_profile(s)
radial_control = [0 0.4 0.8 1 1 1];
height_control = [1 1 1 0.8 0.4 0];
[radial,radial_derivative] = bezier_values(radial_control,s);
[height,height_derivative] = bezier_values(height_control,s);
end


function [value,derivative] = bezier_values(control,s)
degree = numel(control)-1;
value = zeros(size(s));
for k = 0:degree
    basis = nchoosek(degree,k)*(1-s).^(degree-k).*s.^k;
    value = value+control(k+1)*basis;
end
derivative_control = degree*diff(control);
derivative = zeros(size(s));
for k = 0:degree-1
    basis = nchoosek(degree-1,k)*(1-s).^(degree-1-k).*s.^k;
    derivative = derivative+derivative_control(k+1)*basis;
end
end


function [position,derivative] = chunkie_panel_values(outline,panel,q)
t = 2*q-1;
[polynomials,polynomial_derivatives] = lege.pols(t,outline.order-1);
coefficients = outline.position_coefficients(:,:,panel);
position = coefficients*reshape(polynomials,outline.order,[]);
derivative = 2*coefficients*reshape( ...
    polynomial_derivatives,outline.order,[]);
end


function [nodes,faces] = triangulate_core(vertices,center_spacing, ...
        side_spacing,side_start)
x_limits = [min(vertices(1,:)) max(vertices(1,:))];
y_limits = [min(vertices(2,:)) max(vertices(2,:))];

if abs(center_spacing-side_spacing)<= ...
        10*eps(max(center_spacing,side_spacing))
    x_grid = x_limits(1):center_spacing:x_limits(2);
    y_grid = y_limits(1):center_spacing:y_limits(2);
    [grid_x,grid_y] = meshgrid(x_grid,y_grid);
    inside = inpolygon(grid_x(:),grid_y(:), ...
        vertices(1,:),vertices(2,:));
    interior = [grid_x(inside).';grid_y(inside).'];
else
    interior = adaptive_core_points(vertices,center_spacing, ...
        side_spacing,side_start,x_limits,y_limits);
end
interior = discard_near_boundary(interior,vertices, ...
    0.3*min(center_spacing,side_spacing));

number_of_vertices = size(vertices,2);
points = [vertices interior].';
constraints = [(1:number_of_vertices).', ...
    [2:number_of_vertices 1].'];
triangulation_object = delaunayTriangulation(points,constraints);
candidate_faces = triangulation_object.ConnectivityList;
centroids = (points(candidate_faces(:,1),:)+ ...
    points(candidate_faces(:,2),:)+points(candidate_faces(:,3),:))/3;
use = inpolygon(centroids(:,1),centroids(:,2), ...
    vertices(1,:),vertices(2,:));
faces = candidate_faces(use,:).';
nodes = points.';

a = nodes(:,faces(1,:));
b = nodes(:,faces(2,:));
c = nodes(:,faces(3,:));
signed_twice_area = (b(1,:)-a(1,:)).*(c(2,:)-a(2,:))- ...
    (b(2,:)-a(2,:)).*(c(1,:)-a(1,:));
nondegenerate = abs(signed_twice_area)>1e-12;
faces = faces(:,nondegenerate);
signed_twice_area = signed_twice_area(nondegenerate);
flip = signed_twice_area<0;
faces([2 3],flip) = faces([3 2],flip);
end


function interior = adaptive_core_points(vertices,center_spacing, ...
        side_spacing,side_start,x_limits,y_limits)
% Use narrow staggered columns near the two x-extremes and wider columns
% through the middle. The same local spacing is used vertically, so the
% resulting Delaunay triangles remain approximately isotropic.
x_scale = max(abs(vertices(1,:)));
candidate_points = zeros(2,0);
column = 0;
x = x_limits(1)+0.5*local_cap_spacing( ...
    x_limits(1),x_scale,center_spacing,side_spacing,side_start);

while x<x_limits(2)
    local_spacing = local_cap_spacing(x,x_scale,center_spacing, ...
        side_spacing,side_start);
    y_offset = 0.5*local_spacing*(1+mod(column,2));
    y = y_limits(1)+y_offset:local_spacing:y_limits(2);
    candidate_points = [candidate_points, ...
        [x*ones(1,numel(y));y]];
    x = x+sqrt(3)/2*local_spacing;
    column = column+1;
end

inside = inpolygon(candidate_points(1,:),candidate_points(2,:), ...
    vertices(1,:),vertices(2,:));
interior = candidate_points(:,inside);
end


function spacing = local_cap_spacing(x,x_scale,center_spacing, ...
        side_spacing,side_start)
side_coordinate = (abs(x)/x_scale-side_start)/(1-side_start);
side_coordinate = min(max(side_coordinate,0),1);
smooth_transition = side_coordinate.^2.*(3-2*side_coordinate);
spacing = center_spacing+(side_spacing-center_spacing).* ...
    smooth_transition;
end


function kept = discard_near_boundary(points,vertices,distance_limit)
keep = true(1,size(points,2));
number_of_vertices = size(vertices,2);
for edge = 1:number_of_vertices
    next_edge = mod(edge,number_of_vertices)+1;
    a = vertices(:,edge);
    edge_vector = vertices(:,next_edge)-a;
    fraction = sum((points-a).*edge_vector,1)/sum(edge_vector.^2);
    fraction = min(max(fraction,0),1);
    projection = a+edge_vector.*fraction;
    keep = keep & vecnorm(points-projection,2,1)>distance_limit;
end
kept = points(:,keep);
end


function values = planar_triangle_values(vertices,height,u,v)
number_of_nodes = numel(u);
xy = vertices(:,1)+(vertices(:,2)-vertices(:,1)).*u+ ...
    (vertices(:,3)-vertices(:,1)).*v;
du_xy = vertices(:,2)-vertices(:,1);
dv_xy = vertices(:,3)-vertices(:,1);
position = [xy;height*ones(1,number_of_nodes)];
derivative_u = repmat([du_xy;0],1,number_of_nodes);
derivative_v = repmat([dv_xy;0],1,number_of_nodes);
normal = normalized_cross(derivative_u,derivative_v);
values = [position;derivative_u;derivative_v;normal];
end


function normal = normalized_cross(derivative_u,derivative_v)
normal = cross(derivative_u,derivative_v,1);
normal = normal./vecnorm(normal,2,1);
end


function error_value = shared_edge_position_error(S,parts,norder)
[error_value,~] = shared_edge_errors(S,parts,norder);
end


function error_value = shared_edge_normal_error(S,parts,norder)
[~,error_value] = shared_edge_errors(S,parts,norder);
end


function [position_error,normal_error] = shared_edge_errors(S,parts,norder)
sample = linspace(-1,1,max(33,4*norder+1));
position_error = 0;
normal_error = 0;
groups = {parts.top_collar,parts.bottom_collar};
edge_coordinates = [1 2];
for profile_panel = 1:parts.number_of_wall_profile_panels
    groups{end+1} = parts.upper_wall_by_profile(profile_panel,:);
    edge_coordinates(end+1) = 2;
    groups{end+1} = parts.lower_wall_by_profile(profile_panel,:);
    edge_coordinates(end+1) = 2;
end
for group_id = 1:numel(groups)
    group = groups{group_id};
    if edge_coordinates(group_id)==1
        current_basis = polytens.lege.pols(norder, ...
            [ones(size(sample));sample]);
        next_basis = polytens.lege.pols(norder, ...
            [-ones(size(sample));sample]);
    else
        current_basis = polytens.lege.pols(norder, ...
            [sample;ones(size(sample))]);
        next_basis = polytens.lege.pols(norder, ...
            [sample;-ones(size(sample))]);
    end
    for panel = 1:numel(group)
        next_panel = mod(panel,numel(group))+1;
        current_values = S.srccoefs{group(panel)}*current_basis;
        next_values = S.srccoefs{group(next_panel)}*next_basis;
        [position_error,normal_error] = update_edge_errors( ...
            current_values,next_values,position_error,normal_error);
    end
end

wall_profile_end = polytens.lege.pols(norder, ...
    [ones(size(sample));sample]);
wall_profile_start = polytens.lege.pols(norder, ...
    [-ones(size(sample));sample]);
for panel = 1:parts.number_of_outline_panels
    for profile_panel = 1:parts.number_of_wall_profile_panels-1
        upper_current = S.srccoefs{parts.upper_wall_by_profile( ...
            profile_panel,panel)}*wall_profile_end;
        upper_next = S.srccoefs{parts.upper_wall_by_profile( ...
            profile_panel+1,panel)}*wall_profile_start;
        [position_error,normal_error] = update_edge_errors( ...
            upper_current,upper_next,position_error,normal_error);

        lower_current = S.srccoefs{parts.lower_wall_by_profile( ...
            profile_panel,panel)}*wall_profile_end;
        lower_next = S.srccoefs{parts.lower_wall_by_profile( ...
            profile_panel+1,panel)}*wall_profile_start;
        [position_error,normal_error] = update_edge_errors( ...
            lower_current,lower_next,position_error,normal_error);
    end
end

upper_wall_cap = polytens.lege.pols(norder, ...
    [-ones(size(sample));sample]);
upper_wall_equator = polytens.lege.pols(norder, ...
    [ones(size(sample));sample]);
lower_wall_equator = polytens.lege.pols(norder, ...
    [-ones(size(sample));sample]);
lower_wall_cap = polytens.lege.pols(norder, ...
    [ones(size(sample));sample]);
top_collar_outer = polytens.lege.pols(norder, ...
    [sample;-ones(size(sample))]);
bottom_collar_outer = polytens.lege.pols(norder, ...
    [-ones(size(sample));sample]);
for panel = 1:parts.number_of_outline_panels
    top_values = S.srccoefs{parts.top_collar(panel)}* ...
        top_collar_outer;
    upper_cap_values = S.srccoefs{ ...
        parts.upper_wall_by_profile(1,panel)}*upper_wall_cap;
    [position_error,normal_error] = update_edge_errors( ...
        top_values,upper_cap_values,position_error,normal_error);

    upper_values = S.srccoefs{parts.upper_wall_by_profile( ...
        end,panel)}*upper_wall_equator;
    lower_values = S.srccoefs{parts.lower_wall_by_profile( ...
        1,panel)}*lower_wall_equator;
    [position_error,normal_error] = update_edge_errors( ...
        upper_values,lower_values,position_error,normal_error);

    lower_cap_values = S.srccoefs{parts.lower_wall_by_profile( ...
        end,panel)}*lower_wall_cap;
    bottom_values = S.srccoefs{parts.bottom_collar(panel)}* ...
        bottom_collar_outer;
    [position_error,normal_error] = update_edge_errors( ...
        lower_cap_values,bottom_values,position_error,normal_error);
end
end


function [position_error,normal_error] = update_edge_errors( ...
        current_values,next_values,position_error,normal_error)
position_error = max(position_error,max(vecnorm( ...
    current_values(1:3,:)-next_values(1:3,:),2,1)));
current_normal = normalized_cross( ...
    current_values(4:6,:),current_values(7:9,:));
next_normal = normalized_cross( ...
    next_values(4:6,:),next_values(7:9,:));
cosine = sum(current_normal.*next_normal,1);
angle = acosd(min(max(cosine,-1),1));
normal_error = max(normal_error,max(angle));
end
