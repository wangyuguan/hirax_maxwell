% Explicit NRCCIE system for two original closed leaves, without local tips.
% A = [A11 A12; A21 A22], target leaf first, source leaf second.
% Each node carries interleaved [ju;jv;rho]; only A11/A22 contain jump*I.
% Low-order discretization for inspecting the matrix, not a converged solve.
% Matrix and near corrections are built in target-row batches. The full A
% is stored on disk, not allocated in RAM. Read a block from a full run with:
%   M = matfile(output_file);
%   A12 = M.A(leaf_dof_ids{1},leaf_dof_ids{2});
% Set NRCCIE_TWO_LEAF_ROWS_ONLY=1 to build/verify only sampled target rows;
% that mode writes a separate *_rows.mat with their global matrix_rows.
% For full assembly followed by SVD on a server, run run_two_leaf_system_svd.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Parameters: original shape, coarse panelization and order-1 quadrature
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

%% Initialize the output; matrix_rows maps saved rows to the global system
output_dir = fullfile(root,'data');
if ~isfolder(output_dir), mkdir(output_dir); end
suffix = '';
if rows_only, suffix = '_rows'; end
output_file = fullfile(output_dir,sprintf('two_leaf_system_%d_%d_order%d%s.mat', ...
    settings.leaf_ids(1),settings.leaf_ids(2),settings.norder,suffix));
assembly = struct('complete',false,'rows_only',rows_only,'rows_written',0, ...
    'matrix_size',[numel(matrix_rows),ndof],'full_matrix_size',[ndof,ndof], ...
    'quadrature_seconds',0,'matrix_seconds',0,'cross_correction_nonzeros',0);
validation = struct('passed',false,'relative_error',0,'disk_error',0);
save(output_file,'S','parts','settings','leaf_node_ids','leaf_dof_ids', ...
    'matrix_rows','assembly','validation','-v7.3')
M = matfile(output_file,'Writable',true);
% Allocate the complex dataset on disk. All entries are overwritten below.
M.A(numel(matrix_rows),ndof) = complex(NaN,NaN);
probe = randn(ndof,1)+1i*randn(ndof,1);
probe = probe/norm(probe);
opts = struct('zk',settings.zk,'alpha',settings.alpha,'jump',settings.jump);
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = settings.eps_fmm;

%% Stream the corrected self/cross blocks of the combined surface
for first = 1:settings.target_batch_nodes:numel(node_order)
    last = min(first+settings.target_batch_nodes-1,numel(node_order));
    nodes = node_order(first:last);
    rows = reshape(3*nodes+(-2:0).',1,[]);
    saved_rows = 3*first-2:3*last;
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
    M.A(saved_rows,1:ndof) = Ablock;
    assembly.matrix_seconds = assembly.matrix_seconds+toc(timer);
    if any(ismember(nodes,sample_nodes))
        reference = nrccie_apply(S,probe,C,fmm_opts);
        % Only the current target rows have corrections in C.
        err = norm(Ablock*probe-reference(rows))/norm(reference(rows));
        validation.relative_error = max(validation.relative_error,err);
        disk_rows = M.A(saved_rows(1:3),1:ndof);
        err = norm(disk_rows-Ablock(1:3,:),'fro')/norm(Ablock(1:3,:),'fro');
        validation.disk_error = max(validation.disk_error,err);
    end
    assembly.rows_written = saved_rows(end);
    M.assembly = assembly;
    fprintf('Saved matrix rows %d / %d.\n',saved_rows(end),numel(matrix_rows))
    clear Ablock C Crows
end
assembly.complete = true;
validation.passed = validation.relative_error<1e-7 && validation.disk_error<1e-13;
M.assembly = assembly;
M.validation = validation;
fprintf('Dense rows/FMM %.3e; disk round-trip %.3e; cross correction nnz %d.\n', ...
    validation.relative_error,validation.disk_error,assembly.cross_correction_nonzeros)
fprintf('Saved %s\n',output_file)
assert(validation.passed,'System matrix verification failed.')
disp('TWO_LEAF_SYSTEM_MATRIX_PASS')
