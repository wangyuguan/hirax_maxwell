
clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

leaf_radius = 69.25;
thickness = 0.01*leaf_radius;
surface_order = 6;
chunkie_order = 20;
chunkie_n0 = 3;
chunkie_nchs = 3;
chunkie_newton_iterations = 30;
rim_width = 0.028128271246*leaf_radius;
cap_collar_width = rim_width;
cap_mesh_spacing_center = 0.2*leaf_radius;
cap_mesh_spacing_side = cap_mesh_spacing_center;
cap_mesh_side_start = 0.62;
outline_refinement = 2;
wall_profile_refinement = 3;

geometry_options = struct();
geometry_options.norder = surface_order;
geometry_options.chunkie_order = chunkie_order;
geometry_options.chunkie_n0 = chunkie_n0;
geometry_options.chunkie_nchs = chunkie_nchs;
geometry_options.chunkie_newton_iterations = chunkie_newton_iterations;
geometry_options.rim_width = rim_width;
geometry_options.cap_collar_width = cap_collar_width;
geometry_options.cap_mesh_spacing_center = cap_mesh_spacing_center;
geometry_options.cap_mesh_spacing_side = cap_mesh_spacing_side;
geometry_options.cap_mesh_side_start = cap_mesh_side_start;
geometry_options.outline_refinement = outline_refinement;
geometry_options.wall_profile_refinement = wall_profile_refinement;

[S,parts] = hirax_chunkie_leaf_plate_surfer( ...
    thickness,geometry_options);

%% Chunkie smoothing

outline_samples_per_panel = 101;
outline_parameter = linspace(0,1,outline_samples_per_panel);
outline_polynomials = lege.pols( ...
    2*outline_parameter-1,parts.chunkie_order-1);
smooth_outline = zeros(2,parts.number_of_outline_panels* ...
    (outline_samples_per_panel-1)+1);
outline_cursor = 1;
for panel = 1:parts.number_of_outline_panels
    outline_coefficients = ...
        parts.outline_position_coefficients(:,:,panel);
    outline_values = outline_coefficients*reshape( ...
        outline_polynomials,parts.chunkie_order,[]);
    outline_ids = outline_cursor:outline_cursor+ ...
        outline_samples_per_panel-2;
    smooth_outline(:,outline_ids) = outline_values(:,1:end-1);
    outline_cursor = outline_cursor+outline_samples_per_panel-1;
end
smooth_outline(:,end) = smooth_outline(:,1);
raw_polygon = parts.chunkie_raw_vertices(:,[1:end 1]);

figure_handle = figure;
mesh_axes = axes(figure_handle);

outline_display_scale = 1/parts.inner_scale;
S_display = affine_transf(S,diag([ ...
    outline_display_scale outline_display_scale 1]));
plot_surfer_patch_boundaries(mesh_axes,S_display,[0 0 0],0.35, ...
    parts.top_collar,[0;0;0],9);

outline_z = 0.5*thickness;
hold(mesh_axes,'on')
plot3(mesh_axes,smooth_outline(1,:),smooth_outline(2,:), ...
    outline_z*ones(1,size(smooth_outline,2)), ...
    'r-','LineWidth',1.4)
plot3(mesh_axes,raw_polygon(1,:),raw_polygon(2,:), ...
    outline_z*ones(1,size(raw_polygon,2)), ...
    'k--','LineWidth',1.0)
polygon_connecting_points = [ ...
    parts.design.bottom_left,parts.design.bottom_right, ...
    parts.design.right_circle_start,parts.design.right_ellipse_join, ...
    parts.design.left_ellipse_join,parts.design.left_circle_end];
scatter3(mesh_axes,polygon_connecting_points(1,:), ...
    polygon_connecting_points(2,:), ...
    outline_z*ones(1,size(polygon_connecting_points,2)), ...
    28,'k','filled')
hold(mesh_axes,'off')
axis(mesh_axes,'equal')
axis(mesh_axes,'tight')
view(mesh_axes,0,90)
grid(mesh_axes,'off')
zoom(figure_handle,'on')
