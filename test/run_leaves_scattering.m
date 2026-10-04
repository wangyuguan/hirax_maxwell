% Scattering matrix of the complete four-identical-leaf HIRAX dish.
% Scat = Pout * A^(-1) * Pin includes all interactions between the leaves.
% Run on a server: matlab -batch "run('test/run_leaves_scattering.m')"
% For geometry only, set NRCCIE_LEAVES_SCATTERING_GEOMETRY_ONLY=1 first.
% Reuse the saved result for plane waves without solving the BIE again:
%   b = nrccie_planewave_rhs(Sin,settings.zk,settings.alpha,rotations);
%   qout = Scat*b;
%   fields = nrccie_matrix(Sout,targets,settings.zk)*qout;
% fields contains [Ex;Ey;Ez;Hx;Hy;Hz] per target outside the proxy surface.
% qout is an equivalent density, not the physical density restricted to Jout.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Parameters: same original leaves and physical convention as run_leaves_nrccie
settings.norder = 4;
settings.leaf_radius = 69.25; % mm
settings.thickness = 0.6925;
settings.geometry_options = leaf_geometry_options(struct('norder',settings.norder));
settings.geometry_revision = leaf_geometry_revision();
settings.zk = 2*pi/(10*settings.leaf_radius);
settings.alpha = 1;
settings.jump = 0.5;
settings.eps_quad = 1e-11;
settings.eps_fmm = 1e-10;
settings.factor_tol = 1e-8;
settings.skeleton_tol = 1e-7;
settings.max_leaf_nodes = 100; % FLAM box occupancy, not physical leaf size
settings.max_patches = 16; % connected patches in each skeletonization part
settings.quad_batch_size = 128;
settings.rhs_batch_size = 16;
settings.nrotations = 256; % two transverse polarizations per SO(3) sample
settings.outgoing_proxy_scale = 3;
settings.outgoing_proxy_order = 8;
settings.validation_nrotations = 2;
settings.validation_ndirections = 24;
settings.gmres_restart = 50;
settings.gmres_maxit = 5;
settings.gmres_tol = 1e-9;
settings.random_seed = 89;
rng(settings.random_seed)
geometry_only = strcmp(getenv('NRCCIE_LEAVES_SCATTERING_GEOMETRY_ONLY'),'1');
output_dir = fullfile(root,'data');
if ~isfolder(output_dir), mkdir(output_dir); end
output_file = fullfile(output_dir,sprintf('four_leaf_scattering_order%d.mat',settings.norder));

%% Four identical leaves: TL/TR/BL/BR, with full top, bottom and wall surfaces
timer = tic;
[S,parts] = generate_leaves_surfer(settings.thickness,settings.geometry_options);
settings.geometry_options = parts.options;
leaves = parts.leaves;
reflection = leaf_midplane_permutation(parts.base_surfer,parts.base_parts);
geometry = struct('base',parts.base_surfer,'rotation',parts.rotation,'shift',parts.shift);
leaf_node_ids = parts.node_ids;
leaf_patch_ids = parts.patch_ids;
parts = rmfield(parts,{'base_surfer','leaves'}); % keep one compact geometry copy
timings.geometry = toc(timer);
ndof = 3*S.npts;
assert(all(isfinite(S.r),'all') && all(S.wts>0))
fprintf('Four identical leaves, order %d: %d nodes / leaf, %d total, %d DOFs.\n', ...
    settings.norder,geometry.base.npts,S.npts,ndof)
fprintf('Wavelength %.4f mm; hypothetical dense A %.2f GiB (not assembled).\n', ...
    2*pi/settings.zk,16*ndof^2/2^30)
if geometry_only
    geometry_file = fullfile(output_dir,sprintf( ...
        'four_leaf_scattering_geometry_order%d.mat',settings.norder));
    save(geometry_file,'S','geometry','parts','settings','timings', ...
        'leaf_node_ids','leaf_patch_ids','-v7.3')
    fprintf('Geometry-only check saved %s\n',geometry_file)
    return
end

%% Global incoming/outgoing skeletons, built before the large FLAM factor
rotations = sample_so3(settings.nrotations);
fprintf('Incoming skeleton: %d SO(3) samples, %d RHS columns.\n', ...
    settings.nrotations,2*settings.nrotations)
timer = tic;
[Jin,blocks_in,Pmerge_in,info_in] = hierarchical_nrccie_incoming_skeleton( ...
    S,settings.zk,settings.alpha,rotations,settings.skeleton_tol,settings.max_patches);
timings.incoming = toc(timer);
proxy_center = (min(S.r,[],2)+max(S.r,[],2))/2;
dish_extent = max(vecnorm(S.r-proxy_center));
proxy_radius = settings.outgoing_proxy_scale*dish_extent;
proxy = geometries.sphere(proxy_radius,1,proxy_center,settings.outgoing_proxy_order,11);
assert(min(vecnorm(proxy.r-proxy_center))>dish_extent)
fprintf('Outgoing skeleton: %d proxy nodes enclosing the entire dish.\n',proxy.npts)
timer = tic;
[Jout,blocks_out,Pmerge_out,info_out] = hierarchical_nrccie_outgoing_skeleton( ...
    S,proxy,settings.zk,settings.skeleton_tol,settings.max_patches);
timings.outgoing = toc(timer);
incoming_nodes = ceil(Jin(1:3:end)/3);
outgoing_nodes = ceil(Jout(1:3:end)/3);
Sin = node_subset(S,incoming_nodes);
Sout = node_subset(S,outgoing_nodes);
skeleton_counts = zeros(4,2);
for leaf = 1:4
    skeleton_counts(leaf,:) = [nnz(ismember(incoming_nodes,leaf_node_ids{leaf})), ...
        nnz(ismember(outgoing_nodes,leaf_node_ids{leaf}))];
    fprintf('Leaf %d: %d incoming / %d outgoing skeleton nodes.\n', ...
        leaf,skeleton_counts(leaf,1),skeleton_counts(leaf,2))
end
fprintf('Total skeletons: %d incoming / %d outgoing nodes, %d connected parts.\n', ...
    Sin.npts,Sout.npts,numel(blocks_in))

%% Near quadrature: reuse C4 rotations and z reflection of complete leaves
timer = tic;
[corrections,correction_info] = fourfold_quad_corr_mats( ...
    leaves{1},leaves{2},leaves{3},leaves{4},settings.eps_quad,settings.zk, ...
    settings.quad_batch_size,reflection);
timings.quadrature = toc(timer);
storage = whos('corrections');
correction_info.stored_bytes = storage.bytes;
clear leaves
% FLAM block queries/proxy support need four indexable sparse matrices.
% This expands ONLY the near corrections, not the dense BIE system.
timer = tic;
C = cell(1,4);
[C{:}] = fourfold_expand_quad_corr(corrections);
timings.expand_corrections = toc(timer);
storage = whos('C');
correction_info.expanded_bytes = storage.bytes;
fprintf('Near correction storage: compact %.2f GiB, expanded sparse %.2f GiB.\n', ...
    correction_info.stored_bytes/2^30,correction_info.expanded_bytes/2^30)

%% Factor the coupled four-leaf system; FLAM samples matrix blocks on demand
opts = struct('zk',settings.zk,'alpha',settings.alpha,'jump',settings.jump, ...
    'fmm',false,'proxy',true,'verb',1);
timer = tic;
F = nrccie_rskelf(S,C,opts,settings.max_leaf_nodes,settings.factor_tol);
timings.factor = toc(timer);
storage = whos('F');
factor_bytes = storage.bytes;
clear C % the remaining FMM applications use compact symmetry corrections

%% Scattering map, a few incoming skeleton right-hand sides at a time
timer = tic;
Scat = nrccie_scattering_matrix( ...
    F,blocks_in,Pmerge_in,blocks_out,Pmerge_out,settings.rhs_batch_size);
timings.scattering = toc(timer);
assert(isequal(size(Scat),[3*Sout.npts,3*Sin.npts]))
validation = struct('passed',false);
save(output_file,'Scat','Sin','Sout','Jin','Jout','S','geometry','parts','settings', ...
    'leaf_node_ids','leaf_patch_ids','skeleton_counts', ...
    'blocks_in','Pmerge_in','blocks_out','Pmerge_out','info_in','info_out', ...
    'rotations','proxy','proxy_center','proxy_radius','dish_extent', ...
    'correction_info','factor_bytes','timings','validation','-v7.3')
fprintf('Scat: %d x %d; saved %s\n',size(Scat,1),size(Scat,2),output_file)

%% Held-out waves: compare with the corrected FMM operator on the same mesh
% This validates compression, not discretization error or an analytic PEC
% solution. GMRES uses FLAM as a preconditioner; true FMM residuals are checked.
validation.rotations = sample_so3(settings.validation_nrotations);
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
clear candidate_rhs reconstructed_rhs
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = settings.eps_fmm;
correction_apply = @(d)fourfold_apply_quad_corr(corrections,d);
Afun = @(d)nrccie_apply(S,d,correction_apply,fmm_opts);
probe = randn(ndof,1)+1i*randn(ndof,1);
validation.factor_multiply_error = relative_error(rskelf_mv(F,probe),Afun(probe));
reference_density = complex(zeros(size(B)));
factor_density = rskelf_sv(F,B);
preconditioner = @(d)rskelf_sv(F,d);
timer = tic;
for wave = 1:size(B,2)
    [reference_density(:,wave),flag,~,iter] = gmres(Afun,B(:,wave), ...
        min(settings.gmres_restart,ndof),settings.gmres_tol,settings.gmres_maxit, ...
        preconditioner,[],factor_density(:,wave));
    validation.gmres_flag(wave) = flag;
    validation.gmres_iter(wave,:) = iter;
    validation.reference_residual(wave) = relative_error( ...
        Afun(reference_density(:,wave)),B(:,wave));
    validation.factor_residual(wave) = relative_error(Afun(factor_density(:,wave)),B(:,wave));
end
timings.validation_solve = toc(timer);
validation.factor_density_error = relative_error(factor_density,reference_density);

%% Independent exterior targets and extrapolated far-field amplitudes
target_rotations = sample_so3(settings.validation_ndirections);
validation.directions = reshape(target_rotations(:,3,:),3,[]);
validation.far_radii = [1e4,2e4]*dish_extent;
radii = [4*proxy_radius/3,validation.far_radii];
assert(radii(1)>proxy_radius)
scattering_fields = complex(zeros(6*settings.validation_ndirections,size(B,2),3));
reference_fields = scattering_fields;
equivalent_density = Scat*B(Jin,:);
for ir = 1:3
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
fprintf('Held-out RHS %.3e; FLAM multiplication %.3e; max true FMM residual %.3e.\n', ...
    validation.incoming_error,validation.factor_multiply_error,max(validation.reference_residual))
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
disp('FOUR_LEAF_SCATTERING_PASS')

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
