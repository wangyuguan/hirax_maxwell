% SERVER RUN: assemble the explicit NRCCIE system for two original closed
% leaves, without local tips, and save its singular values. From the
% repository root:
%   matlab -batch "run('test/run_two_leaf_system_svd.m')"
% On Polar, submit test/run_two_leaf_system_svd_polar.sh instead.
%
% Needs sibling directories ../fmm3dbie-hirax-dev and ../chunkie, with their
% MATLAB paths and server-native MEX dependencies.
% A = [A11 A12; A21 A22], target leaf first, source leaf second.
% Each node carries interleaved [ju;jv;rho]; only A11/A22 contain jump*I.
% Low-order discretization for inspecting the matrix, not a converged solve.
% Defaults use original leaves 1 and 2; set the patch order in settings.norder
% (order 1: 17712 scalar unknowns, order 3: 67488).
%
% A is assembled in RAM in target-row batches and never written to disk
% (4.67 GiB at order 1, 67.87 GiB at order 3); svd(A) needs about one more
% copy as workspace.
% Only the singular values and a spectrum plot are saved, to
% data/*_svd.mat and *_svd.png. These are the raw coefficient-matrix
% singular values in the Euclidean DOF norm; no quadrature-weighted L2
% rescaling is applied.
%
% NRCCIE_TWO_LEAF_ROWS_ONLY=1 is a local smoke test: it builds/verifies only
% sampled target rows and saves their singular values to a separate
% *_rows_svd.mat. Those do not describe the full system. setenv persists for
% the MATLAB session.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Parameters: original shape, coarse panelization and low-order patches
settings.leaf_ids = [1,2]; % TL/TR/BL/BR = 1/2/3/4
settings.leaf_radius = 69.25; % mm
settings.thickness = 0.6925;
settings.surface_gap = 4;
settings.norder = 1; % do not use 0: constant position patches degenerate
settings.zk = 2*pi/(10*settings.leaf_radius);
settings.alpha = 1;
settings.jump = 0.5;
settings.eps_quad = 1e-8;
settings.eps_fmm = 1e-10;
settings.target_batch_nodes = 32;
settings.geometry_options = struct('norder',settings.norder, ...
    'chunkie_order',20,'chunkie_n0',3,'chunkie_nchs',3, ...
    'chunkie_newton_iterations',30, ...
    'rim_width',0.028128271246*settings.leaf_radius, ...
    'cap_collar_width',0.028128271246*settings.leaf_radius, ...
    'cap_mesh_spacing_center',0.5*settings.leaf_radius, ...
    'cap_mesh_spacing_side',0.5*settings.leaf_radius, ...
    'outline_refinement',1,'wall_profile_refinement',1);
assert(numel(settings.leaf_ids)==2 && numel(unique(settings.leaf_ids))==2 ...
    && all(ismember(settings.leaf_ids,1:4)))
assert(settings.norder>=1 && settings.norder==fix(settings.norder))
rows_only = strcmp(getenv('NRCCIE_TWO_LEAF_ROWS_ONLY'),'1');
rng(79)

%% One original leaf shape, placed twice in the dish
[base,parts] = leaf_plate_surfer(settings.thickness,settings.geometry_options);
settings.geometry_options = parts.options;
center = parts.design.leaf_center_raw;
center(2) = center(2)+(settings.surface_gap-2.83*sqrt(2))/sqrt(2);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
leaves = cell(1,2);
leaf_node_ids = cell(1,2);
leaf_dof_ids = cell(1,2);
sample_nodes = [];
sample_patches = [parts.top_core(1),parts.top_collar(1), ...
    parts.upper_wall(1),parts.lower_wall(1),parts.bottom_collar(1),parts.bottom_core(1)];
for leaf = 1:2
    angle = angles(settings.leaf_ids(leaf));
    R = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
    settings.rotation(:,:,leaf) = R;
    settings.shift(:,leaf) = R*[center;0];
    leaves{leaf} = affine_transf(base,R,settings.shift(:,leaf));
    leaf_node_ids{leaf} = (leaf-1)*base.npts+(1:base.npts);
    leaf_dof_ids{leaf} = reshape(3*leaf_node_ids{leaf}+(-2:0).',1,[]);
    sample_nodes = [sample_nodes,(leaf-1)*base.npts+ ...
        reshape(base.ixyzs(sample_patches),1,[])]; %#ok<AGROW>
end
% Also sample a facing pair, to exercise cross-leaf near quadrature.
[~,i1] = min(vecnorm(leaves{1}.r-mean(leaves{2}.r,2)));
[~,i2] = min(vecnorm(leaves{2}.r-leaves{1}.r(:,i1)));
[~,i1] = min(vecnorm(leaves{1}.r-leaves{2}.r(:,i2)));
sample_nodes = unique([sample_nodes,i1,base.npts+i2]);
S = merge([leaves{:}]);
ndof = 3*S.npts;
node_order = 1:S.npts;
if rows_only, node_order = sample_nodes; end
matrix_rows = reshape(3*node_order+(-2:0).',1,[]);
fprintf('Leaves %s, order %d: %d + %d nodes; full A is %d x %d (%.2f GiB).\n', ...
    mat2str(settings.leaf_ids),settings.norder,base.npts,base.npts, ...
    ndof,ndof,16*ndof^2/2^30)

%% Assemble the corrected self/cross blocks of the combined surface in RAM
output_dir = fullfile(root,'data');
if ~isfolder(output_dir), mkdir(output_dir); end
suffix = '';
if rows_only, suffix = '_rows'; end
stem = sprintf('two_leaf_system_%d_%d_order%d%s', ...
    settings.leaf_ids(1),settings.leaf_ids(2),settings.norder,suffix);
assembly = struct('rows_only',rows_only,'matrix_size',[numel(matrix_rows),ndof], ...
    'full_matrix_size',[ndof,ndof],'quadrature_seconds',0,'matrix_seconds',0, ...
    'cross_correction_nonzeros',0);
validation = struct('passed',false,'relative_error',0);
A = zeros(numel(matrix_rows),ndof,'like',1i);
probe = randn(ndof,1)+1i*randn(ndof,1);
probe = probe/norm(probe);
opts = struct('zk',settings.zk,'alpha',settings.alpha,'jump',settings.jump);
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = settings.eps_fmm;
for first = 1:settings.target_batch_nodes:numel(node_order)
    last = min(first+settings.target_batch_nodes-1,numel(node_order));
    nodes = node_order(first:last);
    rows = reshape(3*nodes+(-2:0).',1,[]);
    block_rows = 3*first-2:3*last;
    timer = tic;
    Crows = nrccie_quad_corr_block(S,settings.eps_quad,settings.zk,[], ...
        settings.target_batch_nodes,nodes);
    % Lift only these target rows into the global indexing expected by the
    % self block evaluator. All other correction rows stay empty/sparse.
    scatter_rows = sparse(nodes,1:numel(nodes),1,S.npts,numel(nodes));
    C = cell(1,4);
    for kernel = 1:4
        C{kernel} = scatter_rows*Crows{kernel};
        for target_leaf = 1:2
            target_rows = ismember(nodes,leaf_node_ids{target_leaf});
            assembly.cross_correction_nonzeros = assembly.cross_correction_nonzeros ...
                +nnz(Crows{kernel}(target_rows,leaf_node_ids{3-target_leaf}));
        end
    end
    assembly.quadrature_seconds = assembly.quadrature_seconds+toc(timer);
    timer = tic;
    Ablock = nrccie_system_block(S,[],C,opts,rows,1:ndof);
    assert(all(isfinite(Ablock),'all'),'Nonfinite system matrix entries.')
    A(block_rows,:) = Ablock;
    assembly.matrix_seconds = assembly.matrix_seconds+toc(timer);
    if any(ismember(nodes,sample_nodes))
        reference = nrccie_apply(S,probe,C,fmm_opts);
        % Only the current target rows have corrections in C.
        err = norm(Ablock*probe-reference(rows))/norm(reference(rows));
        validation.relative_error = max(validation.relative_error,err);
    end
    fprintf('Assembled matrix rows %d / %d.\n',block_rows(end),numel(matrix_rows))
    clear Ablock C Crows
end
validation.passed = validation.relative_error<1e-7;
fprintf('Dense rows/FMM %.3e; cross correction nnz %d.\n', ...
    validation.relative_error,assembly.cross_correction_nonzeros)
assert(validation.passed,'System matrix verification failed.')
disp('TWO_LEAF_SYSTEM_MATRIX_PASS')
clear probe reference

%% Singular values only; A is discarded afterwards
matrix_size = size(A);
scope = 'full system';
if rows_only, scope = 'sampled rows only'; end
fprintf('SVD of %s: %d x %d complex matrix.\n',scope,matrix_size)
timer = tic;
singular_values = svd(A);
svd_seconds = toc(timer);
clear A
relative_singular_values = singular_values/singular_values(1);
fprintf('SVD time %.2f s; sigma_max %.8e; sigma_min %.8e.\n', ...
    svd_seconds,singular_values(1),singular_values(end))
if rows_only
    fprintf('These are singular values of selected rows, not the full system.\n')
end
svd_file = fullfile(output_dir,[stem,'_svd.mat']);
save(svd_file,'singular_values','relative_singular_values', ...
    'svd_seconds','matrix_size','matrix_rows','settings','assembly', ...
    'validation','scope')
fig = figure('Visible','off','Color','w');
semilogy(1:numel(singular_values),singular_values,'o-','MarkerSize',3)
xlabel('Index'); ylabel('Singular value'); grid on
title(sprintf('Two-leaf NRCCIE: %s (%d x %d)',scope,matrix_size))
plot_file = fullfile(output_dir,[stem,'_svd.png']);
exportgraphics(fig,plot_file,'Resolution',180)
close(fig)
fprintf('Saved %s\nSaved %s\n',svd_file,plot_file)
