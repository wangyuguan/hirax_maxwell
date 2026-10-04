% Quadrature symmetry on actual cap, collar and wall patches of a leaf.
% The fixture is an open patch subset, not a closed-leaf convergence test.
% At most 840 nodes enter quadrature; the full leaf selects actual near pairs.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(97)

settings.norder = 4; % Match the four-leaf scattering script's default order.
settings.leaf_radius = 69.25;
settings.thickness = 0.6925;
settings.eps_quad = 1e-9;
settings.zk = pi/settings.leaf_radius;
settings.batch_size = 37; % Deliberately cuts through patches.
geometry_options = leaf_geometry_options(struct('norder',settings.norder, ...
    'cap_mesh_spacing',0.5*settings.leaf_radius,'wall_profile_refinement',3));
[full_base,parts] = leaf_plate_surfer(settings.thickness,geometry_options);
full_reflection = leaf_midplane_permutation(full_base,parts);

% Pick nearby wall panels, then use the actual near list to select cap
% patches that require cross-leaf quadrature at the unchanged physical gap.
center = parts.design.leaf_center_raw;
center(2) = center(2)+(geometry_options.surface_gap-2.83*sqrt(2))/sqrt(2);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
rotations = zeros(3,3,4);
for leaf = 1:4
    a = angles(leaf);
    rotations(:,:,leaf) = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
end
outline_centers = full_base.cms(:,parts.upper_wall_by_profile(1,:));
left = rotations(:,:,1)*(outline_centers+[center;0]);
right = rotations(:,:,2)*(outline_centers+[center;0]);
distance_squared = squeeze(sum((reshape(left,3,[],1)- ...
    reshape(right,3,1,[])).^2,1));
[~,closest] = min(distance_squared(:));
[panel_a,panel_b] = ind2sub(size(distance_squared),closest);
collar = parts.top_collar(panel_a);
core = slicesurfer(full_base,parts.top_core);
target_core = affine_transf(core,rotations(:,:,1),rotations(:,:,1)*[center;0]);
source_core = affine_transf(core,rotations(:,:,2),rotations(:,:,2)*[center;0]);
near = getnear(source_core,target_core);
assert(~isempty(near.col_ind),'The physical adjacent leaves have no near cap pair.')
near_targets = repelem((1:target_core.npts).',diff(near.row_ptr));
near_ratio = vecnorm(target_core.r(:,near_targets)- ...
    source_core.cms(:,near.col_ind),2,1).'./source_core.rads(near.col_ind);
[settings.near_ratio,strongest] = min(near_ratio);
core_patches = unique(parts.top_core([target_core.patch_id(near_targets(strongest)), ...
    near.col_ind(strongest)]),'stable');
settings.near_core_patches = core_patches;
assert(settings.near_ratio<near.rfac)
upper_patches = [reshape(core_patches,1,[]),collar, ...
    parts.upper_wall_by_profile(1,panel_a), ...
    parts.upper_wall_by_profile(end,panel_b)];
lower_patches = full_base.patch_id(full_reflection(full_base.ixyzs(upper_patches)));
selected_patches = [upper_patches,reshape(lower_patches,1,[])];
node_lists = cell(1,numel(selected_patches));
for j = 1:numel(selected_patches)
    patch = selected_patches(j);
    node_lists{j} = full_base.ixyzs(patch):full_base.ixyzs(patch+1)-1;
end
original_nodes = [node_lists{:}];
base = slicesurfer(full_base,selected_patches);
inverse_nodes = zeros(1,full_base.npts);
inverse_nodes(original_nodes) = 1:base.npts;
reflection = inverse_nodes(full_reflection(original_nodes));
assert(all(reflection>0) && isequal(reflection(reflection),1:base.npts))
assert(all(reflection~=(1:base.npts)))
assert(base.npts==2*(numel(core_patches)*(settings.norder+1)*(settings.norder+2)/2+ ...
    3*(settings.norder+1)^2))
position_tol = 2e4*eps(max(1,max(abs(base.r),[],'all')));
assert(max(abs(base.r(:,reflection)-[1;1;-1].*base.r),[],'all')<position_tol)
assert(max(abs(base.n(:,reflection)-[1;1;-1].*base.n),[],'all')<2e4*eps)
assert(max(abs(base.wts(reflection)-base.wts))<2e4*eps(max(base.wts)))
leaves = cell(1,4);
for leaf = 1:4
    R = rotations(:,:,leaf);
    leaves{leaf} = affine_transf(base,R,R*[center;0]);
end
S = merge([leaves{:}]);
assert(S.npts<=840)
settings.geometry_options = geometry_options;
settings.selected_patches = selected_patches;
settings.outline_panels = [panel_a,panel_b];
settings.npts = S.npts;
clear full_base full_reflection inverse_nodes core source_core target_core near
fprintf(['Actual-leaf patch fixture: %d nodes; outline panels %d, %d; ' ...
    'near cap patches %s (distance/radius %.3f).\n'], ...
    S.npts,panel_a,panel_b,mat2str(core_patches),settings.near_ratio)

% Independent reference: library quadrature on every merged target, with no
% symmetry reconstruction and no nrccie_quad_corr_block batching wrapper.
direct = cell(1,4);
direct{1} = em3d.slp.get_quad_corr_mat(S,settings.eps_quad,settings.zk);
[direct{2:4}] = em3d.sgrad.get_quad_corr_mat(S,settings.eps_quad,settings.zk);
[rotation_only,rotation_info] = fourfold_quad_corr_mats(leaves{:}, ...
    settings.eps_quad,settings.zk,settings.batch_size);
[combined,combined_info] = fourfold_quad_corr_mats(leaves{:}, ...
    settings.eps_quad,settings.zk,settings.batch_size,reflection);
assert(any(rotation_info.nonempty(2:4)) && any(combined_info.nonempty(2:4)), ...
    'The selected real patches must exercise cross-leaf near corrections.')
expanded = cell(2,4);
[expanded{1,:}] = fourfold_expand_quad_corr(rotation_only);
[expanded{2,:}] = fourfold_expand_quad_corr(combined);

% Check all sixteen target/source blocks. Gradient errors share a vector
% norm denominator, so a vanishing Cartesian component cannot hide a defect
% or generate a spurious infinite relative error.
errors.matrix = zeros(2,4);
errors.blocks = zeros(2,4,4,4);
errors.reflection = zeros(1,4);
global_reflection = reshape(reflection(:)+(0:3)*base.npts,1,[]);
signs = [1,1,1,-1];
for k = 1:4
    reference_norm = norm(direct{k},'fro');
    errors.reflection(k) = norm(direct{k}(global_reflection,global_reflection)- ...
        signs(k)*direct{k},'fro')/max(realmin,reference_norm);
    for method = 1:2
        errors.matrix(method,k) = norm(expanded{method,k}-direct{k},'fro')/ ...
            max(realmin,reference_norm);
    end
end
for target_leaf = 1:4
    rows = (target_leaf-1)*base.npts+(1:base.npts);
    for source_leaf = 1:4
        cols = (source_leaf-1)*base.npts+(1:base.npts);
        block_norms = cellfun(@(C)norm(C(rows,cols),'fro'),direct);
        scale = [block_norms(1),repmat(norm(block_norms(2:4)),1,3)];
        for method = 1:2
            for k = 1:4
                errors.blocks(method,k,target_leaf,source_leaf) = ...
                    norm(expanded{method,k}(rows,cols)-direct{k}(rows,cols),'fro')/ ...
                    max(realmin,scale(k));
            end
        end
    end
end

% Independent, unrelated complex values on every leaf and scalar channel.
density = randn(4,S.npts)+1i*randn(4,S.npts);
applied = cell(2,4);
[applied{1,:}] = fourfold_apply_quad_corr(rotation_only,density);
[applied{2,:}] = fourfold_apply_quad_corr(combined,density);
errors.callback = zeros(2,4);
for method = 1:2
    for k = 1:4
        expected = (direct{k}*density.').';
        errors.callback(method,k) = norm(applied{method,k}-expected,'fro')/ ...
            max(realmin,norm(expected,'fro'));
    end
end
probe = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
callback = @(d)fourfold_apply_quad_corr(combined,d);
opts = struct('zk',settings.zk,'alpha',1,'fmm',true,'eps_fmm',1e-10);
errors.nrccie = zeros(1,2);
for j = 1:2
    opts.jump = (-1)^j*0.5;
    expected = nrccie_apply(S,probe,direct,opts);
    actual = nrccie_apply(S,probe,callback,opts);
    errors.nrccie(j) = norm(actual-expected)/norm(expected);
end
fprintf('Rotation-only correction errors: %s\n',mat2str(errors.matrix(1,:),5))
fprintf('Rotation + reflection errors:    %s\n',mat2str(errors.matrix(2,:),5))
fprintf('Direct reflection parity errors: %s\n',mat2str(errors.reflection,5))
fprintf('Worst per-block %.3e; callback %.3e; NRCCIE %.3e.\n', ...
    max(errors.blocks,[],'all'),max(errors.callback,[],'all'),max(errors.nrccie))
passed = max([errors.matrix(:);errors.blocks(:);errors.reflection(:); ...
    errors.callback(:);errors.nrccie(:)])<100*settings.eps_quad;
output_dir = fullfile(root,'diagnostics','quadrature_symmetry');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(fullfile(output_dir,'leaf_patch_results.mat'),'settings','errors','passed')
assert(passed,'Real-leaf patch quadrature symmetry checks failed; inspect saved errors.')
fprintf('LEAF_PATCH_QUAD_SYMMETRY_PASS\n')
