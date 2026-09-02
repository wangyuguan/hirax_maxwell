% Four identical leaves or unit spheres in a clover arrangement.

clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('../FMM3D/matlab')
addpath('src')

%% Geometry choice

use_sphere = false;
if_solve = false;
leaf_radius = 69.25;
sphere_radius = 1;
norder = 4;
sphere_subdivisions = 4;
sphere_iptype = 1;

if use_sphere
    geometry_radius = sphere_radius;
    geometry_name = 'unit sphere';
else
    geometry_radius = leaf_radius;
    geometry_name = 'HIRAX leaf';
end

thickness = 0.01*leaf_radius;

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

opts = struct();
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

%% NRCCIE and GMRES parameters

wavelength = 10*geometry_radius;
zk = 2*pi/wavelength;
alpha = 1;

eps_quad = 1e-7;
eps_fmm = 1e-7;
eps_gmres = 1e-6;
gmres_restart = 30;
gmres_restart_cycles = 10;

source_info = struct();
source_info.r = [0;0;10*geometry_radius];
source_info.edips = -geometry_radius^3*[1;1i;0];

%% One master surface mesh

if use_sphere
    S0 = geometries.sphere( ...
        sphere_radius,sphere_subdivisions,[0;0;0], ...
        norder,sphere_iptype);
else
    [S0,parts0] = hirax_chunkie_leaf_plate_surfer(thickness,opts);
    outline_xy = sample_outline(parts0,201);
end

%% Four rotations

theta1 = pi/4;
theta2 = -pi/4;
theta3 = 3*pi/4;
theta4 = -3*pi/4;

R1 = [cos(theta1),-sin(theta1),0; ...
      sin(theta1), cos(theta1),0; ...
                0,           0,1];
R2 = [cos(theta2),-sin(theta2),0; ...
      sin(theta2), cos(theta2),0; ...
                0,           0,1];
R3 = [cos(theta3),-sin(theta3),0; ...
      sin(theta3), cos(theta3),0; ...
                0,           0,1];
R4 = [cos(theta4),-sin(theta4),0; ...
      sin(theta4), cos(theta4),0; ...
                0,           0,1];

%% Translate the four inner edges to a gap of R/30

surface_gap = geometry_radius/30;

if use_sphere
    sphere_center_offset = sphere_radius+surface_gap/2;
    shift1 = [-sphere_center_offset; sphere_center_offset;0];
    shift2 = [ sphere_center_offset; sphere_center_offset;0];
    shift3 = [-sphere_center_offset;-sphere_center_offset;0];
    shift4 = [ sphere_center_offset;-sphere_center_offset;0];
else
    outline_xy1 = R1(1:2,1:2)*outline_xy;
    outline_xy2 = R2(1:2,1:2)*outline_xy;
    outline_xy3 = R3(1:2,1:2)*outline_xy;
    outline_xy4 = R4(1:2,1:2)*outline_xy;

    shift1 = [ ...
        -surface_gap/2-max(outline_xy1(1,:)); ...
         surface_gap/2-min(outline_xy1(2,:)); ...
         0];
    shift2 = [ ...
         surface_gap/2-min(outline_xy2(1,:)); ...
         surface_gap/2-min(outline_xy2(2,:)); ...
         0];
    shift3 = [ ...
        -surface_gap/2-max(outline_xy3(1,:)); ...
        -surface_gap/2-max(outline_xy3(2,:)); ...
         0];
    shift4 = [ ...
         surface_gap/2-min(outline_xy4(1,:)); ...
        -surface_gap/2-max(outline_xy4(2,:)); ...
         0];
end

S1 = affine_transf(S0,R1,shift1);
S2 = affine_transf(S0,R2,shift2);
S3 = affine_transf(S0,R3,shift3);
S4 = affine_transf(S0,R4,shift4);
S = merge([S1,S2,S3,S4]);

if use_sphere
    top_gap = shift2(1)-shift1(1)-2*sphere_radius;
    bottom_gap = shift4(1)-shift3(1)-2*sphere_radius;
    left_gap = shift1(2)-shift3(2)-2*sphere_radius;
    right_gap = shift2(2)-shift4(2)-2*sphere_radius;
else
    outline_xy1 = outline_xy1+shift1(1:2);
    outline_xy2 = outline_xy2+shift2(1:2);
    outline_xy3 = outline_xy3+shift3(1:2);
    outline_xy4 = outline_xy4+shift4(1:2);

    top_gap = min(outline_xy2(1,:))-max(outline_xy1(1,:));
    bottom_gap = min(outline_xy4(1,:))-max(outline_xy3(1,:));
    left_gap = min(outline_xy1(2,:))-max(outline_xy3(2,:));
    right_gap = min(outline_xy2(2,:))-max(outline_xy4(2,:));
end

fprintf('Four identical %s meshes\n',geometry_name)
fprintf('  geometry radius: %.8g\n',geometry_radius)
fprintf('  requested gap R/30: %.8g\n',surface_gap)
fprintf('  top / bottom gaps: %.8g / %.8g\n',top_gap,bottom_gap)
fprintf('  left / right gaps: %.8g / %.8g\n',left_gap,right_gap)
fprintf('  patches per object: %d, nodes per object: %d\n', ...
    S0.npatches,S0.npts)
fprintf('  total patches: %d, total nodes: %d\n',S.npatches,S.npts)

%% Mesh display with visual-only thickness exaggeration

if use_sphere
    z_display_scale = 1;
else
    z_display_scale = 30;
end
S_display = affine_transf(S,diag([1 1 z_display_scale]));

figure(1)
clf
mesh_axes = axes;
plot_surfer_patch_boundaries(mesh_axes,S_display,[0 0 0],0.55, ...
    1:S_display.npatches,[0;0;0],17);
axis(mesh_axes,'equal')
axis(mesh_axes,'tight')
view(mesh_axes,35,32)
grid(mesh_axes,'off')
box(mesh_axes,'on')

if if_solve

%% NRCCIE right-hand side

surface_normal = S.n;
surface_ru = S.du./vecnorm(S.du,2,1);
surface_rv = cross(surface_normal,surface_ru,1);

[einc,hinc] = em3d.incoming_sources( ...
    zk,source_info,S,'electric dipole');

normal_einc = sum(surface_normal.*einc,1);
nxhinc = cross(surface_normal,hinc,1);
nxnxeinc = surface_normal.*normal_einc-einc;
tangent_rhs = nxhinc-alpha*nxnxeinc;

rhs_components = complex(zeros(3,S.npts));
rhs_components(1,:) = sum(surface_ru.*tangent_rhs,1);
rhs_components(2,:) = sum(surface_rv.*tangent_rhs,1);
rhs_components(3,:) = normal_einc;
rhs = rhs_components(:);

%% Self and cross-leaf quadrature corrections

quadrature_timer = tic;
Cslp = em3d.slp.get_quad_corr_mat(S,eps_quad,zk);
[Cx,Cy,Cz] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk);
quadrature_time = toc(quadrature_timer);

fprintf('  quadrature corrections: %.2f s\n',quadrature_time)
fprintf('  correction nonzeros S / dx / dy / dz: %d / %d / %d / %d\n', ...
    nnz(Cslp),nnz(Cx),nnz(Cy),nnz(Cz))

component_node_count = S0.npts;
slp_correction_blocks = zeros(4,4);
for target_component = 1:4
    target_nodes = (target_component-1)*component_node_count+ ...
        (1:component_node_count);
    for source_component = 1:4
        source_nodes = ...
            (source_component-1)*component_node_count+ ...
            (1:component_node_count);
        slp_correction_blocks(target_component,source_component) = ...
            nnz(Cslp(target_nodes,source_nodes));
    end
end
fprintf(['  S correction blocks; rows are target surfaces and columns ' ...
    'are source surfaces:\n'])
disp(slp_correction_blocks)

%% Matrix-free NRCCIE solve with ordinary restarted GMRES

operator = struct();
operator.npts = S.npts;
operator.r = S.r;
operator.wts = S.wts(:).';
operator.n = surface_normal;
operator.ru = surface_ru;
operator.rv = surface_rv;
operator.zk = zk;
operator.alpha = alpha;
operator.eps_fmm = eps_fmm;
operator.Cslp = Cslp;
operator.Cx = Cx;
operator.Cy = Cy;
operator.Cz = Cz;

matvec = @(density) apply_nrccie(density,operator);

solve_timer = tic;
[solution,gmres_flag,gmres_relative_residual,gmres_iterations, ...
    gmres_residual_history] = gmres( ...
    matvec,rhs,gmres_restart,eps_gmres,gmres_restart_cycles);
solve_time = toc(solve_timer);

if gmres_iterations(1)==0
    number_of_iterations = gmres_iterations(2);
else
    number_of_iterations = ...
        (gmres_iterations(1)-1)*gmres_restart+gmres_iterations(2);
end

true_relative_residual = norm(matvec(solution)-rhs)/norm(rhs);

density_components = reshape(solution,3,S.npts);
surface_current = surface_ru.*density_components(1,:)+ ...
    surface_rv.*density_components(2,:);
surface_charge = density_components(3,:);

fprintf('  GMRES iterations: %d, flag: %d\n', ...
    number_of_iterations,gmres_flag)
fprintf('  reported / recomputed residual: %.3e / %.3e\n', ...
    gmres_relative_residual,true_relative_residual)
fprintf('  solve time: %.2f s\n',solve_time)

%% Save data

settings = struct();
settings.use_sphere = use_sphere;
settings.if_solve = if_solve;
settings.geometry_name = geometry_name;
settings.geometry_radius = geometry_radius;
settings.leaf_radius = leaf_radius;
settings.sphere_radius = sphere_radius;
settings.sphere_subdivisions = sphere_subdivisions;
settings.sphere_iptype = sphere_iptype;
settings.thickness = thickness;
settings.surface_gap = surface_gap;
settings.norder = norder;
settings.wavelength = wavelength;
settings.zk = zk;
settings.alpha = alpha;
settings.eps_quad = eps_quad;
settings.eps_fmm = eps_fmm;
settings.eps_gmres = eps_gmres;
settings.gmres_restart = gmres_restart;
settings.gmres_restart_cycles = gmres_restart_cycles;
settings.source_info = source_info;
settings.rotations = cat(3,R1,R2,R3,R4);
settings.shifts = [shift1,shift2,shift3,shift4];

solver_data = struct();
solver_data.gmres_flag = gmres_flag;
solver_data.gmres_iterations = gmres_iterations;
solver_data.number_of_iterations = number_of_iterations;
solver_data.gmres_relative_residual = gmres_relative_residual;
solver_data.true_relative_residual = true_relative_residual;
solver_data.gmres_residual_history = gmres_residual_history;
solver_data.quadrature_time = quadrature_time;
solver_data.solve_time = solve_time;
solver_data.correction_nonzeros = [ ...
    nnz(Cslp),nnz(Cx),nnz(Cy),nnz(Cz)];
solver_data.slp_correction_blocks = slp_correction_blocks;

if ~isfolder('data')
    mkdir('data')
end
if use_sphere
    output_file = sprintf( ...
        'data/run_multiple_leaves_sphere_order%d.mat',norder);
else
    output_file = sprintf('data/run_multiple_leaves_order%d.mat',norder);
end
save(output_file,'S','settings','solver_data','einc','hinc', ...
    'rhs_components','surface_current','surface_charge','-v7.3')
fprintf('  saved %s\n',output_file)

end


function xy = sample_outline(parts,nq)

q = linspace(0,1,nq);
pols = lege.pols(2*q-1,parts.chunkie_order-1);
npan = parts.number_of_outline_panels;
xy = zeros(2,npan*(nq-1));

for ipan = 1:npan
    values = parts.outline_position_coefficients(:,:,ipan)* ...
        reshape(pols,parts.chunkie_order,[]);
    ids = (ipan-1)*(nq-1)+(1:nq-1);
    xy(:,ids) = values(:,1:end-1);
end
end


function y = apply_nrccie(x,operator)

npts = operator.npts;
density_components = reshape(x,3,npts);
surface_current = operator.ru.*density_components(1,:)+ ...
    operator.rv.*density_components(2,:);
surface_charge = density_components(3,:);
density = [surface_current;surface_charge];

source = struct();
source.sources = operator.r;
source.nd = 4;
source.charges = density.*operator.wts;

fmm_output = hfmm3d(operator.eps_fmm,operator.zk,source,2);
potential = reshape(fmm_output.pot,4,npts);
gradient = reshape(fmm_output.grad,4,3,npts);

potential = potential+(operator.Cslp*density.').';
gradient_x = (operator.Cx*density.').';
gradient_y = (operator.Cy*density.').';
gradient_z = (operator.Cz*density.').';
gradient(:,1,:) = gradient(:,1,:)+reshape(gradient_x,4,1,npts);
gradient(:,2,:) = gradient(:,2,:)+reshape(gradient_y,4,1,npts);
gradient(:,3,:) = gradient(:,3,:)+reshape(gradient_z,4,1,npts);

slp_current = potential(1:3,:);
slp_charge = potential(4,:);
gradient_slp_charge = reshape(gradient(4,:,:),3,npts);

curl_slp_current = complex(zeros(3,npts));
curl_slp_current(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,npts);
curl_slp_current(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,npts);
curl_slp_current(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,npts);
divergence_slp_current = reshape( ...
    gradient(1,1,:)+gradient(2,2,:)+gradient(3,3,:),1,npts);

electric_field = 1i*operator.zk*slp_current-gradient_slp_charge;
nxh = cross(operator.n,curl_slp_current,1);
normal_electric_field = sum(operator.n.*electric_field,1);
nxnxe = operator.n.*normal_electric_field-electric_field;

principal_value = complex(zeros(3,npts));
tangent_equation = -nxh+operator.alpha*nxnxe;
principal_value(1,:) = sum(operator.ru.*tangent_equation,1);
principal_value(2,:) = sum(operator.rv.*tangent_equation,1);
principal_value(3,:) = -normal_electric_field+operator.alpha*( ...
    divergence_slp_current-1i*operator.zk*slp_charge);

y = 0.5*density_components+principal_value;
y = y(:);
end
