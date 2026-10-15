% Use actual leaf patch pairs, arbitrary complex densities, and full assembly.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
clear pth dir
[geometry,metadata] = build_leaf_geometry(6);
% Match the production quadrature tolerance; keep the acceptance threshold.
rng(15); zk = .009; eps_quad = 1e-11;
[U,g] = sample_piece(geometry.common,metadata.common_groups,2);
up = leaf_midplane_permutation(U,g);
common = cell(1,4); local = cell(1,4); reflection = cell(1,8);
for j = 1:4
    [D,h] = sample_piece(geometry.local_tips{j},metadata.local_groups{j},j);
    reflection{2*j-1} = up;
    reflection{2*j} = leaf_midplane_permutation(D,h);
    common{j} = affine_transf(U,geometry.rotation(:,:,j),geometry.shift(:,j));
    local{j} = affine_transf(D,geometry.rotation(:,:,j),geometry.shift(:,j));
end
pieces = reshape([common;local],1,[]); S = merge([pieces{:}]);
[C,info] = common_local_quad_corr_mats(common,local,eps_quad,zk,37,reflection);
full = nrccie_quad_corr_block(S,eps_quad,zk,[],37);
d = randn(4,S.npts)+1i*randn(4,S.npts);
y = cell(1,4); [y{:}] = common_local_apply_quad_corr(C,d);
errors = zeros(1,4);
for k = 1:4
    ref = (full{k}*d.').';
    errors(k) = norm(y{k}-ref,'fro')/norm(ref,'fro');
end
assert(max(errors)<1e-7,'Reflected quadrature differs from full assembly.')
% The diagnostic expansion must retain the public full-matrix interface.
expanded = cell(1,4); [expanded{:}] = fourfold_expand_quad_corr(C.common);
base = cell(1,4); dc = d(:,C.common_nodes);
[base{:}] = fourfold_apply_quad_corr(C.common,dc);
for k = 1:4
    ref = (expanded{k}*dc.').';
    assert(norm(base{k}-ref,'fro')<1e-11*norm(ref,'fro'))
end
op = struct('npts',S.npts,'r',S.r,'wts',S.wts(:).','n',S.n, ...
    'ru',S.du./vecnorm(S.du,2,1),'zk',zk,'alpha',1,'eps_fmm',1e-11);
op.rv = cross(S.n,op.ru,1);
op.apply_corrections = @(v)common_local_apply_quad_corr(C,v);
refop = op; refop.apply_corrections = @(v)full_apply(full,v);
x = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
ref = apply_nrccie(x,refop);
operator_error = norm(apply_nrccie(x,op)-ref)/norm(ref);
fprintf('Correction errors %.3e %.3e %.3e %.3e; NRCCIE error %.3e\n',errors,operator_error)
geometry_revision = metadata.geometry_revision;
save(fullfile(root,'diagnostics','midplane_quad_corr_operator_check.mat'), ...
    'geometry_revision','eps_quad','errors','operator_error');
assert(operator_error<1e-7,'NRCCIE relative error %.3e exceeds 1e-7.',operator_error)
[C0,info0] = common_local_quad_corr_mats(common,local,eps_quad,zk,37);
plainop = op; plainop.apply_corrections = @(v)common_local_apply_quad_corr(C0,v);
plain = apply_nrccie(x,plainop);
mirror_error = norm(apply_nrccie(x,op)-plain)/norm(plain);
assert(mirror_error<1e-7)
fprintf('Mirror-only NRCCIE difference %.3e; stored nnz reduction %.3fx\n', ...
    mirror_error,info0.stored_nonzeros/info.stored_nonzeros)
fprintf('Reflected correction errors: %.3e %.3e %.3e %.3e\n',errors)
fprintf('Reflected NRCCIE error: %.3e\n',operator_error)
fprintf('Stored / expanded nonzeros: %d / %d\n',info.stored_nonzeros,info.assembled_nonzeros)
geometry_revision = metadata.geometry_revision;
save(fullfile(root,'diagnostics','midplane_quad_corr_check.mat'), ...
    'geometry_revision','eps_quad','errors','operator_error','mirror_error','info','info0')
disp('MIDPLANE_QUAD_CORR_PASS')

function [T,h] = sample_piece(S,g,ncore)
names = {'top_core','top_collar','upper_wall','lower_wall','bottom_collar','bottom_core'};
selection = {g.top_core(1:ncore),g.top_collar(1), ...
    g.upper_wall_by_profile(:,1).',g.lower_wall_by_profile(:,1).', ...
    g.bottom_collar(1),g.bottom_core(1:ncore)};
ids = []; h = struct();
for j = 1:6
    h.(names{j}) = numel(ids)+(1:numel(selection{j}));
    ids = [ids,selection{j}]; %#ok<AGROW>
end
h.upper_wall_by_profile = h.upper_wall(:);
h.lower_wall_by_profile = h.lower_wall(:);
T = oversample(slicesurfer(S,ids),4);
end

function [p,x,y,z] = full_apply(C,d)
p = (C{1}*d.').'; x = (C{2}*d.').';
y = (C{3}*d.').'; z = (C{4}*d.').';
end
