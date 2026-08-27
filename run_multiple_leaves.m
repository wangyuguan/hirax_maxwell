% Two-leaf geometry and incident-field RHS interpolation diagnostic.

clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

% Use the same temporary single-leaf geometry as run_hirax_single_leaf.m.

dish_radius = 69.25;
thickness = 0.01*dish_radius;
norder = 8;

chunkie_order = 20;
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

opts = [];
opts.norder = norder;
opts.chunkie_order = chunkie_order;
opts.chunkie_n0 = chunkie_n0;
opts.chunkie_nchs = chunkie_nchs;
opts.chunkie_newton_iterations = chunkie_newton_iterations;
opts.rim_width = rim_width;
opts.cap_collar_width = cap_collar_width;
opts.cap_mesh_spacing_center = cap_mesh_spacing_center;
opts.cap_mesh_spacing_side = cap_mesh_spacing_side;
opts.cap_mesh_side_start = cap_mesh_side_start;
opts.outline_refinement = outline_refinement;
opts.wall_profile_refinement = wall_profile_refinement;

[S0,parts0] = hirax_chunkie_dish_plate_surfer(thickness,opts);

% Rotate the rounded crowns outward and the narrow ends inward, as in the
% upper pair of leaf_photos/cloverleaf_width_gerber.png.

theta = pi/4;
R1 = [cos(theta),-sin(theta),0; ...
      sin(theta), cos(theta),0; ...
               0,          0,1];
R2 = [cos(theta), sin(theta),0; ...
     -sin(theta), cos(theta),0; ...
               0,          0,1];

outline_xy = sample_outline(parts0,201);
outline_xy1 = R1(1:2,1:2)*outline_xy;
outline_xy2 = R2(1:2,1:2)*outline_xy;

leaf_gap = dish_radius/30;
shift1 = [-leaf_gap/2-max(outline_xy1(1,:));0;0];
shift2 = [ leaf_gap/2-min(outline_xy2(1,:));0;0];

S1 = affine_transf(S0,R1,shift1);
S2 = affine_transf(S0,R2,shift2);
S = merge([S1,S2]);

outline_xy1 = outline_xy1+shift1(1:2);
outline_xy2 = outline_xy2+shift2(1:2);
sampled_gap = min(outline_xy2(1,:))-max(outline_xy1(1,:));
assert(abs(sampled_gap-leaf_gap) < 1e-12*dish_radius, ...
    'The two-leaf outline gap is incorrect.')

% Far electric point source; lambda=2R gives zk=pi/R.

wavelength = 2*dish_radius;
zk = 2*pi/wavelength;
src_info = [];
src_info.r = [0;0;10*dish_radius];
src_info.edips = -dish_radius^3*[1;1i;0];

fprintf('Two-leaf RHS interpolation diagnostic\n')
fprintf('  radius: %.8g, spacing R/30: %.8g\n', ...
    dish_radius,sampled_gap)
fprintf('  wavelength: %.8g = 2R, zk: %.8g\n',wavelength,zk)
fprintf('  source: [%.8g %.8g %.8g]\n',src_info.r)
fprintf('  patches: %d, points: %d, order: %d\n', ...
    S.npatches,S.npts,norder)

% Native n dot Einc data.

[Einc,~] = em3d.incoming_sources( ...
    zk,src_info,S,'electric dipole');
rhs = sum(S.n.*Einc,1);

% Interpolate to a finer surface and independently reevaluate the analytic
% incident field there. No quadrature correction or solve is performed.

norder_int = min(2*norder,20);
S_int = oversample(S,norder_int);
rhs_int = interpolate_data( ...
    S,rhs,S_int.patch_id,S_int.uvs_targ);

[Einc_int,~] = em3d.incoming_sources( ...
    zk,src_info,S_int,'electric dipole');
rhs_ex = sum(S_int.n.*Einc_int,1);

rhs_abs = abs(rhs_ex);
rhs_int_abs = abs(rhs_int);
err_abs = abs(rhs_ex-rhs_int);
rhs_scale = max(abs(rhs_ex));
err_rel = err_abs/rhs_scale;
err_pos = err_rel(err_rel > 0);
if isempty(err_pos)
    err_floor = realmin;
else
    err_floor = max(min(err_pos),realmin);
end
err_log10 = log10(max(err_rel,err_floor));

fprintf('  interpolation order: %d\n',norder_int)
fprintf('  maximum normalized RHS error: %.3e\n',max(err_rel))

% Display the thin wall with the same vertical exaggeration used by
% test_dish_mesh_rhs.m. The distant source is intentionally outside the
% plots so that it does not collapse the leaf-scale axis limits.

wall_display_scale = 30;
A_display = diag([1 1 wall_display_scale]);
S_display = affine_transf(S,A_display);
S_int_display = affine_transf(S_int,A_display);

figure(1)
clf
rhs_axes = axes;
plot(S_int_display,rhs_abs,'EdgeColor','none')
hold(rhs_axes,'on')
plot_surfer_patch_boundaries(rhs_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
hold(rhs_axes,'off')
axis(rhs_axes,'tight')
view(rhs_axes,35,32)
grid(rhs_axes,'off')
box(rhs_axes,'on')
colorbar(rhs_axes)
title(rhs_axes,'Exact $|\mathbf{n}\cdot\mathbf{E}_{\rm inc}|$', ...
    'Interpreter','latex')

figure(2)
clf
rhs_int_axes = axes;
plot(S_int_display,rhs_int_abs,'EdgeColor','none')
hold(rhs_int_axes,'on')
plot_surfer_patch_boundaries(rhs_int_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
hold(rhs_int_axes,'off')
axis(rhs_int_axes,'tight')
view(rhs_int_axes,35,32)
grid(rhs_int_axes,'off')
box(rhs_int_axes,'on')
colorbar(rhs_int_axes)
title(rhs_int_axes,'Interpolated $|\mathbf{n}\cdot\mathbf{E}_{\rm inc}|$', ...
    'Interpreter','latex')

rhs_clim = [0,max([rhs_abs rhs_int_abs])];
clim(rhs_axes,rhs_clim)
clim(rhs_int_axes,rhs_clim)

figure(3)
clf
err_axes = axes;
plot(S_int_display,err_log10,'EdgeColor','none')
hold(err_axes,'on')
plot_surfer_patch_boundaries(err_axes,S_display,[0 0 0],0.35, ...
    1:S_display.npatches,[0;0;0],9);
hold(err_axes,'off')
axis(err_axes,'tight')
view(err_axes,35,32)
grid(err_axes,'off')
box(err_axes,'on')
clim(err_axes,[min(err_log10),max(err_log10)])
colorbar(err_axes)
title(err_axes, ...
    '$\log_{10}(|e_{\rm rhs}|/\max|\mathbf{n}\cdot\mathbf{E}_{\rm inc}|)$', ...
    'Interpreter','latex')

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
title(mesh_axes,'Two-leaf surface mesh')


function xy = sample_outline(parts,nq)

q = linspace(0,1,nq);
pols = lege.pols(2*q-1,parts.chunkie_order-1);
npan = parts.number_of_outline_panels;
xy = zeros(2,npan*(nq-1));

for ipan = 1:npan
    vals = parts.outline_position_coefficients(:,:,ipan)* ...
        reshape(pols,parts.chunkie_order,[]);
    ids = (ipan-1)*(nq-1)+(1:nq-1);
    xy(:,ids) = vals(:,1:end-1);
end
end
