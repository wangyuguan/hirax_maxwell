% Scattering matrix for one ORIGINAL closed leaf of the HIRAX dish.
% No modified tips or common/local geometry decomposition.
% Scat = Pout * A^(-1) * Pin, with the physical PEC jump +1/2.
% Reuse the result without a BIE solve:
%   b_in = nrccie_planewave_rhs(Sin,settings.zk,settings.alpha,rotations);
%   fields = nrccie_matrix(Sout,targets,settings.zk)*(Scat*b_in);
% fields are interleaved [Ex;Ey;Ez;Hx;Hy;Hz] at each target.
% The original mesh is large even at low order. For geometry inspection ONLY:
%   setenv('NRCCIE_LEAF_SCATTERING_GEOMETRY_ONLY','1'); run this script;
%   setenv('NRCCIE_LEAF_SCATTERING_GEOMETRY_ONLY','');
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Physical parameters from run_leaves_nrccie; original mesh from run_single_leaf
settings.leaf_id = 1; % TL/TR/BL/BR = 1/2/3/4; only this leaf is constructed
settings.leaf_radius = 69.25; % mm
settings.thickness = 0.6925;
settings.surface_gap = 4;
settings.zk = 2*pi/(10*settings.leaf_radius);
settings.alpha = 1;
settings.jump = 0.5;
settings.norder = 4;
settings.geometry_options = struct('norder',settings.norder, ...
    'chunkie_order',20,'chunkie_n0',3,'chunkie_nchs',3, ...
    'chunkie_newton_iterations',30, ...
    'rim_width',0.028128271246*settings.leaf_radius, ...
    'cap_collar_width',0.028128271246*settings.leaf_radius, ...
    'cap_mesh_spacing_center',0.2*settings.leaf_radius, ...
    'cap_mesh_spacing_side',0.2*settings.leaf_radius, ...
    'cap_mesh_side_start',0.62,'outline_refinement',1, ...
    'wall_profile_refinement',3);
settings.eps_quad = 1e-11;
settings.eps_fmm = 1e-10;
settings.factor_tol = 1e-8;
settings.skeleton_tol = 1e-7;
settings.max_leaf_nodes = 100; % FLAM leaf occupancy, not total surface size
settings.max_patches = 16; % each skeletonization group is connected
settings.quad_batch_size = 128;
settings.rhs_batch_size = 16;
settings.nrotations = 128; % two transverse polarizations per SO(3) sample
settings.outgoing_proxy_scale = 3;
settings.outgoing_proxy_order = 8;
settings.gmres_restart = 100;
settings.gmres_maxit = 5;
settings.gmres_tol = 1e-9;
settings.random_seed = 73;
rng(settings.random_seed)
geometry_only = strcmp(getenv('NRCCIE_LEAF_SCATTERING_GEOMETRY_ONLY'),'1');
output_dir = fullfile(root,'data');
if ~isfolder(output_dir), mkdir(output_dir); end
output_file = fullfile(output_dir,sprintf( ...
    'single_leaf_scattering_leaf%d_order%d.mat',settings.leaf_id,settings.norder));

%% Closed original leaf: top, bottom, collar and rounded side wall
timer = tic;
[S,parts] = leaf_plate_surfer(settings.thickness,settings.geometry_options);
% Place this original leaf at the corresponding position in the dish.
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
angle = angles(settings.leaf_id);
rotation = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
center = parts.design.leaf_center_raw;
center(2) = center(2)+(settings.surface_gap-2.83*sqrt(2))/sqrt(2);
shift = rotation*[center;0];
S = affine_transf(S,rotation,shift);
settings.rotation = rotation;
settings.shift = shift;
settings.geometry_options = parts.options;
timings.geometry = toc(timer);
ndof = 3*S.npts;
assert(all(isfinite(S.r),'all') && all(S.wts>0))
assert(~isempty(parts.top) && ~isempty(parts.bottom) && ~isempty(parts.wall))
fprintf('Original leaf %d: %d patches, %d nodes, %d DOFs; no local tips.\n', ...
    settings.leaf_id,S.npatches,S.npts,ndof)
fprintf('Thickness %.4f mm, wavelength %.4f mm; hypothetical dense A %.2f GiB.\n', ...
    settings.thickness,2*pi/settings.zk,16*ndof^2/2^30)
if geometry_only
    geometry_file = fullfile(output_dir,sprintf( ...
        'single_leaf_scattering_geometry_leaf%d_order%d.mat', ...
        settings.leaf_id,settings.norder));
    save(geometry_file,'S','parts','settings','timings','-v7.3')
    fprintf('Geometry-only check: saved %s; no quadrature or solve was run.\n',geometry_file)
    return
end

%% Corrected BIE and FLAM factorization, with its local two-sided proxy
opts = struct('zk',settings.zk,'alpha',settings.alpha,'jump',settings.jump, ...
    'fmm',false,'proxy',true,'verb',1);
timer = tic;
C = nrccie_quad_corr_block(S,settings.eps_quad,settings.zk,[], ...
    settings.quad_batch_size);
timings.quadrature = toc(timer);
timer = tic;
F = nrccie_rskelf(S,C,opts,settings.max_leaf_nodes,settings.factor_tol);
timings.factor = toc(timer);
storage = whos('F');
factor_bytes = storage.bytes;

%% Incoming boundary nodes and outgoing equivalent-density nodes
rotations = sample_so3(settings.nrotations);
timer = tic;
[Jin,blocks_in,Pmerge_in,info_in] = hierarchical_nrccie_incoming_skeleton( ...
    S,settings.zk,settings.alpha,rotations,settings.skeleton_tol,settings.max_patches);
timings.incoming = toc(timer);
proxy_center = (min(S.r,[],2)+max(S.r,[],2))/2;
leaf_extent = max(vecnorm(S.r-proxy_center));
proxy_radius = settings.outgoing_proxy_scale*leaf_extent;
proxy = geometries.sphere(proxy_radius,1,proxy_center,settings.outgoing_proxy_order,11);
assert(min(vecnorm(proxy.r-proxy_center))>leaf_extent)
timer = tic;
[Jout,blocks_out,Pmerge_out,info_out] = hierarchical_nrccie_outgoing_skeleton( ...
    S,proxy,settings.zk,settings.skeleton_tol,settings.max_patches);
timings.outgoing = toc(timer);
Sin = node_subset(S,ceil(Jin(1:3:end)/3));
Sout = node_subset(S,ceil(Jout(1:3:end)/3));
fprintf('Skeletons: %d incoming / %d outgoing complete nodes, %d connected parts.\n', ...
    Sin.npts,Sout.npts,numel(blocks_in))

%% Scattering matrix, formed a few right-hand sides at a time
timer = tic;
Scat = nrccie_scattering_matrix( ...
    F,blocks_in,Pmerge_in,blocks_out,Pmerge_out,settings.rhs_batch_size);
timings.scattering = toc(timer);
assert(isequal(size(Scat),[3*Sout.npts,3*Sin.npts]))
validation = struct('passed',false);
save(output_file,'Scat','Sin','Sout','Jin','Jout','S','settings','parts', ...
    'blocks_in','Pmerge_in','blocks_out','Pmerge_out','info_in','info_out', ...
    'rotations','proxy','proxy_center','proxy_radius','leaf_extent', ...
    'factor_bytes','timings','validation','-v7.3')
fprintf('Scat: %d x %d; saved %s\n',size(Scat,1),size(Scat,2),output_file)

%% Held-out plane wave, both polarizations; independent corrected FMM solve
% This checks compression against the SAME discrete leaf, not discretization
% accuracy or an analytic PEC solution for this shape.
validation.rotations = sample_so3(1);
B = nrccie_planewave_rhs(S,settings.zk,settings.alpha,validation.rotations);
candidate_rhs = Pmerge_in*B(Jin,:);
reconstructed_rhs = complex(zeros(size(B)));
for part = 1:numel(blocks_in)
    b = blocks_in(part);
    dofs = reshape(3*b.nodes+(-2:0).',1,[]);
    reconstructed_rhs(dofs,:) = b.P*candidate_rhs(b.candidate_dofs,:);
end
w = sqrt(repelem(S.wts(:),3));
validation.incoming_error = relative_error(w.*reconstructed_rhs,w.*B);
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = settings.eps_fmm;
Afun = @(density)nrccie_apply(S,density,C,fmm_opts);
probe = randn(ndof,1)+1i*randn(ndof,1);
validation.factor_multiply_error = relative_error(rskelf_mv(F,probe),Afun(probe));
reference_density = complex(zeros(size(B)));
factor_density = rskelf_sv(F,B);
timer = tic;
for wave = 1:size(B,2)
    [reference_density(:,wave),flag,~,iter] = gmres(Afun,B(:,wave), ...
        settings.gmres_restart,settings.gmres_tol,settings.gmres_maxit);
    validation.gmres_flag(wave) = flag;
    validation.gmres_iter(wave,:) = iter;
    validation.reference_residual(wave) = relative_error( ...
        Afun(reference_density(:,wave)),B(:,wave));
    validation.factor_residual(wave) = relative_error( ...
        Afun(factor_density(:,wave)),B(:,wave));
end
timings.validation_solve = toc(timer);
validation.factor_density_error = relative_error(factor_density,reference_density);

%% Fields at independent exterior targets and extrapolated far-field amplitudes
target_rotations = sample_so3(24);
validation.directions = reshape(target_rotations(:,3,:),3,[]);
validation.far_radii = [1e4,2e4]*leaf_extent;
radii = [4*leaf_extent,validation.far_radii];
scattering_fields = complex(zeros(6*24,size(B,2),numel(radii)));
reference_fields = scattering_fields;
equivalent_density = Scat*B(Jin,:); % not the physical density restricted to Jout
for ir = 1:numel(radii)
    targets.r = radii(ir)*validation.directions;
    if ir==1, targets.r = targets.r+proxy_center; end
    scattering_fields(:,:,ir) = nrccie_matrix(Sout,targets,settings.zk)*equivalent_density;
    for part = 1:numel(blocks_out)
        nodes = blocks_out(part).nodes;
        dofs = reshape(3*nodes+(-2:0).',1,[]);
        K = nrccie_matrix(node_subset(S,nodes),targets,settings.zk);
        reference_fields(:,:,ir) = reference_fields(:,:,ir)+K*reference_density(dofs,:);
    end
end
scale = validation.far_radii.*exp(-1i*settings.zk*validation.far_radii);
validation.farfield_scattering = 2*scale(2)*scattering_fields(:,:,3) ...
    -scale(1)*scattering_fields(:,:,2);
validation.farfield_reference = 2*scale(2)*reference_fields(:,:,3) ...
    -scale(1)*reference_fields(:,:,2);
validation.exterior_field_error = relative_error(scattering_fields(:,:,1),reference_fields(:,:,1));
validation.farfield_error = relative_error(validation.farfield_scattering,validation.farfield_reference);
save(output_file,'validation','scattering_fields','reference_fields','timings','-append')
fprintf('Held-out RHS %.3e; FLAM multiplication %.3e, max residual %.3e.\n', ...
    validation.incoming_error,validation.factor_multiply_error,max(validation.factor_residual))
fprintf('Scattering vs FMM/GMRES: exterior %.3e, farfield %.3e.\n', ...
    validation.exterior_field_error,validation.farfield_error)
assert(all(validation.gmres_flag==0) && ...
    max(validation.reference_residual)<5*settings.gmres_tol,'Reference solve did not converge.')
assert(validation.incoming_error<10*settings.skeleton_tol)
assert(max([validation.factor_multiply_error,validation.factor_residual])<100*settings.factor_tol)
assert(max([validation.exterior_field_error,validation.farfield_error]) ...
    <100*settings.factor_tol+200*settings.skeleton_tol)
validation.passed = true;
save(output_file,'validation','-append')
disp('SINGLE_LEAF_SCATTERING_PASS')


function T = node_subset(S,nodes)
T.npts = numel(nodes);
T.r = S.r(:,nodes);
T.dru = S.dru(:,nodes);
T.drv = S.drv(:,nodes);
T.n = S.n(:,nodes);
T.wts = S.wts(nodes);
end


function err = relative_error(value,reference)
err = norm(value-reference,'fro')/max(norm(reference,'fro'),realmin);
end
