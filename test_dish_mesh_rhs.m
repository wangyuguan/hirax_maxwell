
clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

dish_radius = 69.25;
thickness = 0.01*dish_radius;
surface_order = 16;
wall_display_scale = 30;

wavelength = 10*dish_radius;
zk = 2*pi/wavelength;
p0 = dish_radius^3*[1;1i;0];
source_gap = dish_radius/30;
source_height = 0;
outline_samples_per_panel = 2001;

chunkie_order = 10;
chunkie_n0 = 3;
chunkie_nchs = 3;
chunkie_newton_iterations = 30;
rim_width = 0.028128271246*dish_radius;
cap_collar_width = rim_width;

cap_mesh_spacing_center = 0.2*dish_radius;
cap_mesh_spacing_side = cap_mesh_spacing_center;
cap_mesh_side_start = 0.62;
outline_refinement = 1;
wall_profile_refinement = 3;

geometry_options = struct();
geometry_options.norder = surface_order;
geometry_options.chunkie_order = chunkie_order;
geometry_options.chunkie_n0 = chunkie_n0;
geometry_options.chunkie_nchs = chunkie_nchs;
geometry_options.chunkie_newton_iterations = ...
    chunkie_newton_iterations;
geometry_options.rim_width = rim_width;
geometry_options.cap_collar_width = cap_collar_width;
geometry_options.cap_mesh_spacing_center = cap_mesh_spacing_center;
geometry_options.cap_mesh_spacing_side = cap_mesh_spacing_side;
geometry_options.cap_mesh_side_start = cap_mesh_side_start;
geometry_options.outline_refinement = outline_refinement;
geometry_options.wall_profile_refinement = wall_profile_refinement;

[S,parts] = hirax_chunkie_dish_plate_surfer( ...
    thickness,geometry_options);

outline_parameter = linspace(0,1,outline_samples_per_panel);
outline_polynomials = lege.pols( ...
    2*outline_parameter-1,parts.chunkie_order-1);
source_anchor = [-inf;NaN;0];
for panel = 1:parts.number_of_outline_panels
    outline_coefficients = ...
        parts.outline_position_coefficients(:,:,panel);
    outline_values = outline_coefficients*reshape( ...
        outline_polynomials,parts.chunkie_order,[]);
    [candidate_x,candidate_index] = max(outline_values(1,:));
    if candidate_x>source_anchor(1)
        source_anchor = [outline_values(:,candidate_index);0];
    end
end
source_point = source_anchor+[source_gap;0;source_height];

%% Native normal electric field

source_info = struct();
source_info.r = source_point;
source_info.edips = -p0;
[einc,~] = em3d.incoming_sources( ...
    zk,source_info,S,'electric dipole');

normal_einc = sum(S.n.*einc,1);

%% Interpolate to a finer surface grid

interpolation_order = min(2*surface_order,20);
S_interpolation = oversample(S,interpolation_order);
normal_einc_interpolated = interpolate_data( ...
    S,normal_einc,S_interpolation.patch_id,S_interpolation.uvs_targ);

[einc_interpolation,~] = em3d.incoming_sources( ...
    zk,source_info,S_interpolation,'electric dipole');

normal_einc_exact = sum(S_interpolation.n.*einc_interpolation,1);
normal_einc_magnitude = abs(normal_einc_exact);
interpolated_normal_einc_magnitude = abs(normal_einc_interpolated);
absolute_difference = abs( ...
    normal_einc_exact-normal_einc_interpolated);
normal_einc_maximum = max(abs(normal_einc_exact));
normalized_absolute_difference = ...
    absolute_difference/normal_einc_maximum;
positive_normalized_difference = ...
    normalized_absolute_difference(normalized_absolute_difference>0);
difference_display_floor = min(positive_normalized_difference);
log10_normalized_absolute_difference = log10(max( ...
    normalized_absolute_difference,difference_display_floor));

%% Display surfaces

display_transform = diag([1 1 wall_display_scale]);
S_display = affine_transf(S,display_transform);
S_interpolation_display = affine_transf( ...
    S_interpolation,display_transform);
source_point_display = display_transform*source_point;

%% Normal electric field

figure(1)
clf
rhs_axes = axes;
plot(S_interpolation_display,normal_einc_magnitude,'EdgeColor','none')
hold(rhs_axes,'on')
plot_surfer_patch_boundaries(rhs_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
scatter3(rhs_axes,source_point_display(1),source_point_display(2), ...
    source_point_display(3),48,'r','filled')
hold(rhs_axes,'off')
axis(rhs_axes,'tight')
view(rhs_axes,35,32)
grid(rhs_axes,'off')
box(rhs_axes,'on')
colorbar(rhs_axes)

%% Interpolated normal electric field

figure(2)
clf
interpolated_rhs_axes = axes;
plot(S_interpolation_display,interpolated_normal_einc_magnitude, ...
    'EdgeColor','none')
hold(interpolated_rhs_axes,'on')
plot_surfer_patch_boundaries( ...
    interpolated_rhs_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
scatter3(interpolated_rhs_axes,source_point_display(1), ...
    source_point_display(2),source_point_display(3),48,'r','filled')
hold(interpolated_rhs_axes,'off')
axis(interpolated_rhs_axes,'tight')
view(interpolated_rhs_axes,35,32)
grid(interpolated_rhs_axes,'off')
box(interpolated_rhs_axes,'on')
colorbar(interpolated_rhs_axes)

rhs_color_maximum = max([ ...
    normal_einc_magnitude interpolated_normal_einc_magnitude]);
clim(rhs_axes,[0 rhs_color_maximum])
clim(interpolated_rhs_axes,[0 rhs_color_maximum])

%% Normalized absolute difference

figure(3)
clf
difference_axes = axes;
plot(S_interpolation_display,log10_normalized_absolute_difference, ...
    'EdgeColor','none')
hold(difference_axes,'on')
plot_surfer_patch_boundaries( ...
    difference_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
scatter3(difference_axes,source_point_display(1), ...
    source_point_display(2),source_point_display(3),48,'r','filled')
hold(difference_axes,'off')
axis(difference_axes,'tight')
view(difference_axes,35,32)
grid(difference_axes,'off')
box(difference_axes,'on')
difference_color_limits = [ ...
    min(log10_normalized_absolute_difference), ...
    max(log10_normalized_absolute_difference)];
clim(difference_axes,difference_color_limits)
colorbar(difference_axes)

%% Mesh

figure(4)
clf
mesh_axes = axes;
plot_surfer_patch_boundaries(mesh_axes,S_display,[0 0 0],0.55, ...
    1:S_display.npatches,[0;0;0],17);
hold(mesh_axes,'on')
scatter3(mesh_axes,source_point_display(1),source_point_display(2), ...
    source_point_display(3),48,'r','filled')
hold(mesh_axes,'off')
axis(mesh_axes,'equal')
axis(mesh_axes,'tight')
view(mesh_axes,35,32)
grid(mesh_axes,'off')
box(mesh_axes,'on')
