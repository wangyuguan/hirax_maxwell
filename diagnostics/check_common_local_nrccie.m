% Compare compressed corrections and NRCCIE against full-surface assembly.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath(genpath('src'))
addpath('../fmm3dbie-hirax-dev/FMM3D/matlab')
clear pth dir
rng(7)
zk = .7; eps_quad = 1e-8;
master = geometries.sphere(1,1,[0;0;0],3,1);
U = slicesurfer(master,1:6);
D = slicesurfer(master,7:12);
common = cell(1,4); local = cell(1,4);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
for j = 1:4
    R = [cos(angles(j)),-sin(angles(j)),0;sin(angles(j)),cos(angles(j)),0;0,0,1];
    shift = R*[0;1.6;0];
    common{j} = affine_transf(U,R,shift);
    % Unequal local node counts, with shared U-D boundaries on each sphere.
    tip = split_patches(D,1:j-1);
    local{j} = affine_transf(tip,R,shift);
end
pieces = reshape([common;local],1,[]);
S = merge([pieces{:}]);
[C,info] = common_local_quad_corr_mats(common,local,eps_quad,zk,37);
full = cell(1,4);
full{1} = em3d.slp.get_quad_corr_mat(S,eps_quad,zk);
[full{2:4}] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk);
d = randn(4,S.npts)+1i*randn(4,S.npts);
y = cell(1,4);
[y{:}] = common_local_apply_quad_corr(C,d);
errors = zeros(1,4);
for k = 1:4
    ref = (full{k}*d.').';
    errors(k) = norm(y{k}-ref,'fro')/norm(ref,'fro');
end
assert(max(errors)<2e-6,'Partitioned quadrature differs from full assembly.')
op = struct('npts',S.npts,'r',S.r,'wts',S.wts(:).', ...
    'n',S.n,'ru',S.du./vecnorm(S.du,2,1),'zk',zk,'alpha',1,'eps_fmm',1e-10);
op.rv = cross(S.n,op.ru,1);
op.apply_corrections = @(v)common_local_apply_quad_corr(C,v);
refop = op;
refop.apply_corrections = @(v)full_apply(full,v);
x = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
ref = apply_nrccie(x,refop);
operator_error = norm(apply_nrccie(x,op)-ref)/norm(ref);
assert(operator_error<2e-6)
[sol,flag,~,~,resvec] = gmres(@(v)apply_nrccie(v,op),ref,30,1e-8,10);
full_residual = norm(apply_nrccie(sol,refop)-ref)/norm(ref);
assert(flag==0 && full_residual<2e-6)
fprintf('Correction errors: %.3e %.3e %.3e %.3e\n',errors)
fprintf('NRCCIE error %.3e; full residual %.3e; iterations %d\n', ...
    operator_error,full_residual,numel(resvec)-1)
fprintf('Stored / expanded nonzeros: %d / %d\n',info.stored_nonzeros,info.assembled_nonzeros)
disp('COMMON_LOCAL_NRCCIE_PASS')

function [p,x,y,z] = full_apply(C,d)
p = (C{1}*d.').'; x = (C{2}*d.').';
y = (C{3}*d.').'; z = (C{4}*d.').';
end
