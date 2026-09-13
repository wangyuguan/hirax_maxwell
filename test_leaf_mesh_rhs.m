
clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

leaf_radius = 69.25;
thickness = 0.01*leaf_radius;
surface_order = 12;
wall_display_scale = 30;

wavelength = 2*leaf_radius;
zk = 2*pi/wavelength;
p0 = leaf_radius^3*[1;1i;0];

chunkie_order = 20;
chunkie_n0 = 3;
chunkie_nchs = 3;
chunkie_newton_iterations = 30;
rim_width = 0.028128271246*leaf_radius;
cap_collar_width = rim_width;

cap_mesh_spacing_center = 0.2*leaf_radius;
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

[S,parts] = hirax_chunkie_leaf_plate_surfer( ...
    thickness,geometry_options);

% Overhead dipole of run_multiple_leaves.m. The four-leaf right-hand side
% is the incident field restricted to the surface and carries no
% scattering, so every leaf sees this same data up to a scalar phase and a
% single leaf answers the four-leaf question.
source_point = [0;0;10*leaf_radius];

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

%% Per-patch spectral tail of the same function

% surf_fun_error returns the absolute infinity norm of the high-degree
% coefficient tail on every patch, so the input is normalized by its
% global maximum to make the result relative. A per-patch denominator
% would blow up wherever the normal field passes through zero.
normal_einc_scale = max(abs(normal_einc));
patch_error = abs(surf_fun_error(S,normal_einc/normal_einc_scale)).';
positive_patch_error = patch_error(patch_error>0);
patch_error_display_floor = min(positive_patch_error);
log10_patch_error = log10(max(patch_error,patch_error_display_floor));

%% Display surfaces

display_transform = diag([1 1 wall_display_scale]);
S_display = affine_transf(S,display_transform);
S_interpolation_display = affine_transf( ...
    S_interpolation,display_transform);

% The figures deliberately carry no source marker. The source sits at
% z = 10*leaf_radius, which the thickness exaggeration pushes out to
% z = 300*leaf_radius while the leaf spans a few tens of display units.
% Marking it would stretch every axis tight z range by three orders of
% magnitude and flatten the leaf to a line.

%% Normal electric field

figure(1)
clf
rhs_axes = axes;
plot(S_interpolation_display,normal_einc_magnitude,'EdgeColor','none')
hold(rhs_axes,'on')
plot_surfer_patch_boundaries(rhs_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
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
axis(mesh_axes,'equal')
axis(mesh_axes,'tight')
view(mesh_axes,35,32)
grid(mesh_axes,'off')
box(mesh_axes,'on')

%% Per-patch spectral tail

figure(5)
clf
patch_error_axes = axes;
plot(S_display,log10_patch_error,'EdgeColor','none')
hold(patch_error_axes,'on')
plot_surfer_patch_boundaries( ...
    patch_error_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
hold(patch_error_axes,'off')
axis(patch_error_axes,'tight')
view(patch_error_axes,35,32)
grid(patch_error_axes,'off')
box(patch_error_axes,'on')
clim(patch_error_axes, ...
    [min(log10_patch_error) max(log10_patch_error)])
colorbar(patch_error_axes)

%% Summary

weighted_l2_relative_difference = sqrt(sum( ...
    absolute_difference.^2.*S_interpolation.wts(:).'))/sqrt(sum( ...
    abs(normal_einc_exact).^2.*S_interpolation.wts(:).'));

fprintf('HIRAX leaf right-hand-side resolution check\n')
fprintf('  surface order: %d, interpolation order: %d\n', ...
    surface_order,interpolation_order)
fprintf('  patches: %d, nodes: %d, oversampled nodes: %d\n', ...
    S.npatches,S.npts,S_interpolation.npts)
fprintf('  source point: (%.8g, %.8g, %.8g)\n',source_point)
fprintf('  nodal interpolation error, max / weighted L2: %.3e / %.3e\n', ...
    max(normalized_absolute_difference),weighted_l2_relative_difference)
fprintf('  per-patch tail, min / median / max: %.3e / %.3e / %.3e\n', ...
    min(patch_error),median(patch_error),max(patch_error))
fprintf('  patches with tail above 1e-8 / 1e-6: %d / %d\n', ...
    sum(patch_error>1e-8),sum(patch_error>1e-6))
