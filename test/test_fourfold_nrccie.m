% Small regression for rotation/reflection reuse in the NRCCIE callback.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(83)

% Exactly planar boxes keep positions and derivatives consistent at order 1.
% A curved order-1 fixture has unresolved self-PV quadrature unrelated to
% the rotation/reflection indexing tested here. Merge order: TL/TR/BL/BR.
uv = polytens.lege.nodes(1);
half_width = [0.45,0.35,0.2];
values = zeros(12,24);
face = 0;
for normal_axis = 1:3
    tangent_axes = mod(normal_axis+(0:1),3)+1;
    for side = [-1,1]
        face = face+1;
        center = zeros(3,1); center(normal_axis) = side*half_width(normal_axis);
        du = zeros(3,1); du(tangent_axes(1)) = half_width(tangent_axes(1));
        dv = zeros(3,1); dv(tangent_axes(2)) = side*half_width(tangent_axes(2));
        normal = zeros(3,1); normal(normal_axis) = side;
        values(:,4*face-3:4*face) = ...
            [center+du*uv(1,:)+dv*uv(2,:);repmat([du;dv;normal],1,4)];
    end
end
base = surfer(6,1,values,11);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
leaves = cell(1,4);
for j = 1:4
    a = angles(j);
    R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    leaves{j} = affine_transf(base,R,R*[0;1;0]);
end
S = merge([leaves{:}]);
assert(S.npts<=200)

% Match each node with its partner under z -> -z.
reflection = zeros(1,base.npts);
for j = 1:base.npts
    point = [1;1;-1].*base.r(:,j);
    [distance,reflection(j)] = min(vecnorm(base.r-point,2,1));
    assert(distance<1e-12)
end
assert(isequal(reflection(reflection),1:base.npts) && ...
    all(reflection~=(1:base.npts)))

eps_quad = 1e-8;
opts = struct('zk',0.9,'alpha',1,'jump',0.5, ...
    'fmm',true,'eps_fmm',1e-10);
% The reference uses neither symmetry reuse nor our target-batching wrapper.
direct = cell(1,4);
direct{1} = em3d.slp.get_quad_corr_mat(S,eps_quad,opts.zk);
[direct{2:4}] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,opts.zk);
batched = nrccie_quad_corr_block(S,eps_quad,opts.zk,[],7);
batch_error = zeros(1,4);
for k = 1:4
    batch_error(k) = norm(batched{k}-direct{k},'fro')/norm(direct{k},'fro');
end
assert(max(batch_error)<100*eps_quad)

% Independent scalar densities also exercise every Cartesian current channel.
density = randn(4,S.npts)+1i*randn(4,S.npts);
% Unrelated complex [ju,jv,rho] values on every node of all four surfaces.
probe = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
permutations = {[],reflection};
mode_names = {'rotation','rotation + reflection'};
matrix_error = zeros(2,4);
block_error = zeros(4,4,4,2);
apply_error = zeros(2,4);
direct_apply_error = zeros(2,4);
nrccie_error = zeros(2,2);
for mode = 1:2
    [compressed,info] = fourfold_quad_corr_mats( ...
        leaves{1},leaves{2},leaves{3},leaves{4}, ...
        eps_quad,opts.zk,7,permutations{mode});
    expanded = cell(1,4);
    [expanded{:}] = fourfold_expand_quad_corr(compressed);
    % Ensure cross-surface corrections are tested as well as self blocks.
    assert(any(info.nonempty(2:4)))
    applied = cell(1,4);
    [applied{:}] = fourfold_apply_quad_corr(compressed,density);
    for k = 1:4
        matrix_error(mode,k) = norm(expanded{k}-direct{k},'fro')/norm(direct{k},'fro');
        expected = (expanded{k}*density.').';
        apply_error(mode,k) = norm(applied{k}-expected,'fro')/norm(expected,'fro');
        expected = (direct{k}*density.').';
        direct_apply_error(mode,k) = norm(applied{k}-expected,'fro')/norm(expected,'fro');
    end
    for it = 1:4
        rows = (it-1)*base.npts+(1:base.npts);
        for is = 1:4
            cols = (is-1)*base.npts+(1:base.npts);
            % Joint gradient norm avoids division by an identically zero
            % Cartesian component without hiding a weak cross-leaf block.
            gradient_scale = sqrt(sum(cellfun( ...
                @(C)norm(C(rows,cols),'fro')^2,direct(2:4))));
            for k = 1:4
                scale = gradient_scale;
                if k==1, scale = norm(direct{k}(rows,cols),'fro'); end
                difference = norm(expanded{k}(rows,cols)-direct{k}(rows,cols),'fro');
                if scale==0
                    assert(difference==0,'An empty direct block became nonzero.')
                else
                    block_error(it,is,k,mode) = difference/scale;
                end
            end
        end
    end
    callback = @(d)fourfold_apply_quad_corr(compressed,d);
    for side = 1:2
        opts.jump = (-1)^side*0.5;
        via_callback = nrccie_apply(S,probe,callback,opts);
        via_direct = nrccie_apply(S,probe,direct,opts);
        nrccie_error(mode,side) = norm(via_callback-via_direct)/norm(via_direct);
    end
    fprintf('%s: kernels %s; worst block %.3e; NRCCIE %.3e.\n', ...
        mode_names{mode},mat2str(matrix_error(mode,:),5), ...
        max(block_error(:,:,:,mode),[],'all'),max(nrccie_error(mode,:)))
end
assert(max(matrix_error,[],'all')<100*eps_quad)
assert(max(block_error,[],'all')<100*eps_quad)
assert(max(apply_error,[],'all')<1e-12)
assert(max(direct_apply_error,[],'all')<100*eps_quad)
assert(max(nrccie_error,[],'all')<100*eps_quad)
output_dir = fullfile(root,'diagnostics','quadrature_symmetry');
if ~isfolder(output_dir), mkdir(output_dir); end
save(fullfile(output_dir,'four_box_results.mat'),'matrix_error','block_error', ...
    'apply_error','direct_apply_error','nrccie_error','batch_error','eps_quad')
fprintf('FOURFOLD_NRCCIE_PASS\n')
