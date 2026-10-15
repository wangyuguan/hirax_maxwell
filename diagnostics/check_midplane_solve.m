% Closed asymmetric conductors: reflected solve versus independently assembled kernels.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath(genpath('src'))
addpath('../fmm3dbie-hirax-dev/FMM3D/matlab')
clear pth dir
order = 3; zk = .7; eps_quad = 1e-9; eps_gmres = 1e-8;
common = cell(1,4); local = cell(1,4); reflection = cell(1,8);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
U = cube(order,.30);
for j = 1:4
    R = [cos(angles(j)),-sin(angles(j)),0;sin(angles(j)),cos(angles(j)),0;0,0,1];
    common{j} = affine_transf(U,R,R*[0;1.2;0]);
    local{j} = affine_transf(cube(order,.125+.015*j),R,R*[0;.6;0]);
    reflection{2*j-1} = midplane_nodes(common{j});
    reflection{2*j} = midplane_nodes(local{j});
end
pieces = reshape([common;local],1,[]); S = merge([pieces{:}]);
[C,info] = common_local_quad_corr_mats(common,local,eps_quad,zk,2000,reflection);
full = nrccie_quad_corr_block(S,eps_quad,zk,[],2000);
op = struct('npts',S.npts,'r',S.r,'wts',S.wts(:).','n',S.n, ...
    'ru',S.du./vecnorm(S.du,2,1),'zk',zk,'alpha',1,'eps_fmm',1e-11);
op.rv = cross(S.n,op.ru,1);
op.apply_corrections = @(d)common_local_apply_quad_corr(C,d);
refop = op; refop.apply_corrections = @(d)full_apply(full,d);
source_info = struct('r',[3;2;4],'edips',[1;.3i;.2]);
[einc,hinc] = em3d.incoming_sources(zk,source_info,S,'electric dipole');
normal_einc = sum(S.n.*einc,1);
tangent_rhs = cross(S.n,hinc,1)-op.alpha*(S.n.*normal_einc-einc);
rhs = [sum(op.ru.*tangent_rhs,1);sum(op.rv.*tangent_rhs,1);normal_einc];
rhs = rhs(:);
rng(16); x = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
ref = apply_nrccie(x,refop);
operator_error = norm(apply_nrccie(x,op)-ref)/norm(ref);
[solution,flag,relres,iter,resvec] = gmres(@(v)apply_nrccie(v,op),rhs,60,eps_gmres,10);
mirror_residual = norm(apply_nrccie(solution,op)-rhs)/norm(rhs);
full_residual = norm(apply_nrccie(solution,refop)-rhs)/norm(rhs);
fprintf('Closed-cube N %d: operator difference %.3e; GMRES flag %d, iterations %d\n', ...
    S.npts,operator_error,flag,numel(resvec)-1)
fprintf('Mirrored / independently assembled residual: %.3e / %.3e\n', ...
    mirror_residual,full_residual)
save(fullfile(root,'diagnostics','midplane_solve_check.mat'), ...
    'order','zk','eps_quad','eps_gmres','operator_error','flag','relres','iter', ...
    'resvec','mirror_residual','full_residual','info')
assert(operator_error<1e-8 && flag==0 && ...
    mirror_residual<2*eps_gmres && full_residual<2*eps_gmres)
disp('MIDPLANE_SOLVE_PASS')

function S = cube(order,halfwidth)
uv = polytens.lege.nodes(order); np = size(uv,2);
normals = [eye(3),-eye(3)]; vals = zeros(12,6*np);
for f = 1:6
    n = normals(:,f); u = [1;0;0];
    if abs(n(1))==1, u = [0;1;0]; end
    v = cross(n,u); ii = (f-1)*np+(1:np);
    vals(:,ii) = [halfwidth*(n+u*uv(1,:)+v*uv(2,:)); ...
        repmat(halfwidth*u,1,np);repmat(halfwidth*v,1,np);repmat(n,1,np)];
end
S = surfer(6,order,vals,11);
end

function p = midplane_nodes(S)
ref = [1;1;-1].*S.r;
distance = (S.r(1,:).'-ref(1,:)).^2+ ...
    (S.r(2,:).'-ref(2,:)).^2+(S.r(3,:).'-ref(3,:)).^2;
[error,p] = min(distance,[],1);
assert(max(error)<1e-24 && isequal(p(p),1:S.npts) && all(p~=(1:S.npts)))
assert(max(abs(S.n(:,p)-[1;1;-1].*S.n),[],'all')<1e-12)
w = S.wts(:).'; assert(max(abs(w(p)-w))<1e-12)
end

function [p,x,y,z] = full_apply(C,d)
p = (C{1}*d.').'; x = (C{2}*d.').';
y = (C{3}*d.').'; z = (C{4}*d.').';
end
