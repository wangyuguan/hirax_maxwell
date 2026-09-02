% Compare two recursive-skeletonization strategies for the two-sphere NRCCIE.

clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../FLAM/startup.m')
addpath('../FMM3D/matlab')

rng(1)

%% Small two-sphere geometry

a = 1;
d = a/30;
na = [1;2;2];
norder = 3;
iptype = 11;

c = a+d/2;
c0 = [-c,c;0,0;0,0];

abc = a*ones(3,1);
S1 = geometries.ellipsoid(abc,na,c0(:,1),norder,iptype);
S2 = geometries.ellipsoid(abc,na,c0(:,2),norder,iptype);
S = merge([S1,S2]);

n = S.n;
ru = S.du./vecnorm(S.du,2,1);
rv = cross(n,ru,1);

wavelength = 2*a;
zk = 2*pi/wavelength;
alpha = 1;
eps_quad = 1e-10;
eps_fmm = 1e-12;

fprintf('Two-sphere NRCCIE skeletonization test\n')
fprintf('  radius: %.8g, gap: %.8g\n',a,d)
fprintf('  surface nodes: %d, NRCCIE unknowns: %d\n',S.npts,3*S.npts)
fprintf('  patches: %d, surface order: %d\n',S.npatches,norder)

%% Add-subtract quadrature corrections

tquad = tic;
C = cell(1,4);
C{1} = em3d.slp.get_quad_corr_mat(S,eps_quad,zk);
[C{2},C{3},C{4}] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk);
quadrature_time = toc(tquad);

i1 = 1:S1.npts;
i2 = S1.npts+(1:S2.npts);
cross_corrections = 0;
for ikernel = 1:4
    cross_corrections = cross_corrections+nnz(C{ikernel}(i1,i2));
    cross_corrections = cross_corrections+nnz(C{ikernel}(i2,i1));
end

fprintf('  quadrature-correction time: %.3f s\n',quadrature_time)
fprintf('  cross-surface correction entries: %d\n',cross_corrections)

%% NRCCIE matrix generator

geom = [];
geom.r = S.r;
geom.n = n;
geom.ru = ru;
geom.rv = rv;
geom.wts = S.wts(:).';

op = [];
op.npts = S.npts;
op.geom = geom;
op.zk = zk;
op.alpha = alpha;
op.C = C;
op.eps_fmm = eps_fmm;

ndof = 3*S.npts;
ids = 1:ndof;
Afun = @(i,j) nrccie_matrix(i,j,op);

tmatrix = tic;
A = Afun(ids,ids);
matrix_time = toc(tmatrix);
matrix_info = whos('A');

fprintf('  dense reference matrix time: %.3f s\n',matrix_time)
fprintf('  dense reference matrix memory: %.2f MB\n',matrix_info.bytes/2^20)

%% Verify the generated entries against the FMM NRCCIE matvec

v = randn(ndof,1)+1i*randn(ndof,1);
y_dense = A*v;
y_fmm = nrccie_fmm_matvec(v,op);
entry_check = norm(y_dense-y_fmm)/norm(y_fmm);

fprintf('  dense-matrix versus FMM matvec: %.3e\n',entry_check)

%% Incident-field right-hand side and dense reference solution

src_info = [];
src_info.r = a*[0.30;-0.40;2.80];
src_info.edips = [1+0.30i;-0.40+0.20i;0.25-0.55i];

[Einc,Hinc] = em3d.incoming_sources(zk,src_info,S,'electric dipole');
nEinc = sum(n.*Einc,1);
nxHinc = cross(n,Hinc,1);
nxnxEinc = n.*nEinc-Einc;
rhs_t = nxHinc-alpha*nxnxEinc;

rhs = complex(zeros(3,S.npts));
rhs(1,:) = sum(ru.*rhs_t,1);
rhs(2,:) = sum(rv.*rhs_t,1);
rhs(3,:) = nEinc;
rhs = rhs(:);

tdense = tic;
solution_dense = A\rhs;
dense_solve_time = toc(tdense);
dense_residual = norm(A*solution_dense-rhs)/norm(rhs);

fprintf('  dense solve time: %.3f s\n',dense_solve_time)
fprintf('  dense relative residual: %.3e\n',dense_residual)

%% Shared rskelf settings

coordinates = repelem(S.r,1,3);
occupancy = 24;
skeleton_tolerance = 1e-5;

rskelf_options = [];
rskelf_options.symm = 'n';
rskelf_options.Tmax = 2;
rskelf_options.verb = 1;

%% Method 1: direct global off-diagonal compression

fprintf('\nMethod 1: rskelf without a proxy\n')
tfactor1 = tic;
F1 = rskelf(Afun,coordinates,occupancy,skeleton_tolerance,[],rskelf_options);
factor_time1 = toc(tfactor1);
factor_info1 = whos('F1');

rskelf_mv(F1,v);
tmv1 = tic;
y1 = rskelf_mv(F1,v);
mv_time1 = toc(tmv1);

rskelf_sv(F1,rhs);
tsolve1 = tic;
solution1 = rskelf_sv(F1,rhs);
solve_time1 = toc(tsolve1);

mv_error1 = norm(y1-y_dense)/norm(y_dense);
residual1 = norm(A*solution1-rhs)/norm(rhs);
solution_error1 = norm(solution1-solution_dense)/norm(solution_dense);

clear F1 y1 solution1

%% Method 2: one proxy space shared by S, dSdx, dSdy, and dSdz

proxy_count = 96;
proxy_reference = fibonacci_sphere(proxy_count,1.5);
proxyfun = @(x,slf,nbr,l,ctr) joint_kernel_proxy( ...
    x,slf,nbr,l,ctr,proxy_reference,op);

fprintf('\nMethod 2: rskelf with a joint primitive-kernel proxy\n')
tfactor2 = tic;
F2 = rskelf(Afun,coordinates,occupancy,skeleton_tolerance, ...
    proxyfun,rskelf_options);
factor_time2 = toc(tfactor2);
factor_info2 = whos('F2');

rskelf_mv(F2,v);
tmv2 = tic;
y2 = rskelf_mv(F2,v);
mv_time2 = toc(tmv2);

rskelf_sv(F2,rhs);
tsolve2 = tic;
solution2 = rskelf_sv(F2,rhs);
solve_time2 = toc(tsolve2);

mv_error2 = norm(y2-y_dense)/norm(y_dense);
residual2 = norm(A*solution2-rhs)/norm(rhs);
solution_error2 = norm(solution2-solution_dense)/norm(solution_dense);

%% Comparison

method = ["global off-diagonal";"joint kernel proxy"];
factor_seconds = [factor_time1;factor_time2];
factor_megabytes = [factor_info1.bytes;factor_info2.bytes]/2^20;
apply_seconds = [mv_time1;mv_time2];
matrix_apply_error = [mv_error1;mv_error2];
solve_seconds = [solve_time1;solve_time2];
relative_residual = [residual1;residual2];
relative_solution_error = [solution_error1;solution_error2];

results = table(method,factor_seconds,factor_megabytes,apply_seconds, ...
    matrix_apply_error,solve_seconds,relative_residual, ...
    relative_solution_error);

fprintf('\nComparison\n')
disp(results)


function A = nrccie_matrix(i,j,op)

if isempty(i) || isempty(j)
    A = complex(zeros(length(i),length(j)));
    return
end

parts = nrccie_matrix_parts(i,j,op.geom,op.geom, ...
    op.zk,op.alpha,op.C);
A = parts{1}+parts{2}+parts{3}+parts{4};
A = A+0.5*(i(:) == j(:).');
end


function parts = nrccie_matrix_parts(i,j,target,source,zk,alpha,C)

i = i(:);
j = j(:).';

target_node = ceil(i/3);
source_node = ceil(j/3);
target_component = mod(i-1,3)+1;
source_component = mod(j-1,3)+1;

dx = reshape(target.r(1,target_node),[],1)- ...
    reshape(source.r(1,source_node),1,[]);
dy = reshape(target.r(2,target_node),[],1)- ...
    reshape(source.r(2,source_node),1,[]);
dz = reshape(target.r(3,target_node),[],1)- ...
    reshape(source.r(3,source_node),1,[]);
rr = sqrt(dx.^2+dy.^2+dz.^2);

G = complex(zeros(size(rr)));
Gx = complex(zeros(size(rr)));
Gy = complex(zeros(size(rr)));
Gz = complex(zeros(size(rr)));

use = rr > 0;
expikr = exp(1i*zk*rr(use));
G(use) = expikr./(4*pi*rr(use));
radial = (1i*zk*rr(use)-1).*expikr./(4*pi*rr(use).^3);
Gx(use) = radial.*dx(use);
Gy(use) = radial.*dy(use);
Gz(use) = radial.*dz(use);

source_weights = reshape(source.wts(source_node),1,[]);
K = cell(1,4);
K{1} = G.*source_weights;
K{2} = Gx.*source_weights;
K{3} = Gy.*source_weights;
K{4} = Gz.*source_weights;

if ~isempty(C)
    for ikernel = 1:4
        K{ikernel} = K{ikernel}+ ...
            full(C{ikernel}(target_node,source_node));
    end
end

is_ju = source_component == 1;
is_jv = source_component == 2;
is_rho = source_component == 3;

jx = reshape(source.ru(1,source_node),1,[]).*is_ju+ ...
    reshape(source.rv(1,source_node),1,[]).*is_jv;
jy = reshape(source.ru(2,source_node),1,[]).*is_ju+ ...
    reshape(source.rv(2,source_node),1,[]).*is_jv;
jz = reshape(source.ru(3,source_node),1,[]).*is_ju+ ...
    reshape(source.rv(3,source_node),1,[]).*is_jv;
rho = is_rho;

n1 = reshape(target.n(1,target_node),[],1);
n2 = reshape(target.n(2,target_node),[],1);
n3 = reshape(target.n(3,target_node),[],1);
ru1 = reshape(target.ru(1,target_node),[],1);
ru2 = reshape(target.ru(2,target_node),[],1);
ru3 = reshape(target.ru(3,target_node),[],1);
rv1 = reshape(target.rv(1,target_node),[],1);
rv2 = reshape(target.rv(2,target_node),[],1);
rv3 = reshape(target.rv(3,target_node),[],1);

parts = cell(1,4);
for ikernel = 1:4
    Q = K{ikernel};
    Ex = complex(zeros(size(Q)));
    Ey = complex(zeros(size(Q)));
    Ez = complex(zeros(size(Q)));
    curlx = complex(zeros(size(Q)));
    curly = complex(zeros(size(Q)));
    curlz = complex(zeros(size(Q)));
    divj = complex(zeros(size(Q)));
    slp_rho = complex(zeros(size(Q)));

    if ikernel == 1
        Ex = 1i*zk*Q.*jx;
        Ey = 1i*zk*Q.*jy;
        Ez = 1i*zk*Q.*jz;
        slp_rho = Q.*rho;
    elseif ikernel == 2
        Ex = -Q.*rho;
        curly = -Q.*jz;
        curlz = Q.*jy;
        divj = Q.*jx;
    elseif ikernel == 3
        Ey = -Q.*rho;
        curlx = Q.*jz;
        curlz = -Q.*jx;
        divj = Q.*jy;
    else
        Ez = -Q.*rho;
        curlx = -Q.*jy;
        curly = Q.*jx;
        divj = Q.*jz;
    end

    nE = n1.*Ex+n2.*Ey+n3.*Ez;
    nxHx = n2.*curlz-n3.*curly;
    nxHy = n3.*curlx-n1.*curlz;
    nxHz = n1.*curly-n2.*curlx;

    eqx = -nxHx+alpha*(n1.*nE-Ex);
    eqy = -nxHy+alpha*(n2.*nE-Ey);
    eqz = -nxHz+alpha*(n3.*nE-Ez);
    eqrho = -nE+alpha*(divj-1i*zk*slp_rho);

    parts{ikernel} = ...
        (target_component == 1).*(ru1.*eqx+ru2.*eqy+ru3.*eqz)+ ...
        (target_component == 2).*(rv1.*eqx+rv2.*eqy+rv3.*eqz)+ ...
        (target_component == 3).*eqrho;
end
end


function [Kpxy,nbr] = joint_kernel_proxy( ...
    x,slf,nbr,l,ctr,proxy_reference,op)

l = l(:);
ctr = ctr(:);
proxy_points = proxy_reference.*l+ctr;

proxy_normals = (proxy_points-ctr)./(l.^2);
proxy_normals = proxy_normals./vecnorm(proxy_normals,2,1);

frame_seed = repmat([0;0;1],1,size(proxy_points,2));
use_y_axis = abs(proxy_normals(3,:)) > 0.9;
frame_seed(:,use_y_axis) = repmat([0;1;0],1,nnz(use_y_axis));
proxy_ru = cross(frame_seed,proxy_normals,1);
proxy_ru = proxy_ru./vecnorm(proxy_ru,2,1);
proxy_rv = cross(proxy_normals,proxy_ru,1);

proxy_geom = [];
proxy_geom.r = proxy_points;
proxy_geom.n = proxy_normals;
proxy_geom.ru = proxy_ru;
proxy_geom.rv = proxy_rv;
proxy_geom.wts = mean(op.geom.wts)*ones(1,size(proxy_points,2));

proxy_ids = 1:3*size(proxy_points,2);
outgoing = nrccie_matrix_parts(proxy_ids,slf, ...
    proxy_geom,op.geom,op.zk,op.alpha,[]);
incoming = nrccie_matrix_parts(slf,proxy_ids, ...
    op.geom,proxy_geom,op.zk,op.alpha,[]);

Kpxy = complex(zeros(0,length(slf)));
for ikernel = 1:4
    Kpair = [outgoing{ikernel};incoming{ikernel}.'];
    Kscale = norm(Kpair,'fro');
    if Kscale > 0
        Kpair = Kpair/Kscale;
    end
    Kpxy = [Kpxy;Kpair];
end

scaled_distance = (x(:,nbr)-ctr)./l;
nbr = nbr(sum(scaled_distance.^2,1) < 1.5^2);
end


function points = fibonacci_sphere(npoints,radius)

indices = 0:npoints-1;
z = 1-2*(indices+0.5)/npoints;
phi = pi*(3-sqrt(5))*indices;
rxy = sqrt(1-z.^2);
points = radius*[rxy.*cos(phi);rxy.*sin(phi);z];
end


function y = nrccie_fmm_matvec(x,op)

npts = op.npts;
dens_uv = reshape(x,3,npts);
zjvec = op.geom.ru.*dens_uv(1,:)+op.geom.rv.*dens_uv(2,:);
rho = dens_uv(3,:);
dens = [zjvec;rho];

srcinfo = [];
srcinfo.sources = op.geom.r;
srcinfo.nd = 4;
srcinfo.charges = dens.*op.geom.wts;

U = hfmm3d(op.eps_fmm,op.zk,srcinfo,2);
pot = reshape(U.pot,4,npts);
grad = reshape(U.grad,4,3,npts);

pot = pot+(op.C{1}*dens.').';
gx = (op.C{2}*dens.').';
gy = (op.C{3}*dens.').';
gz = (op.C{4}*dens.').';
grad(:,1,:) = grad(:,1,:)+reshape(gx,4,1,npts);
grad(:,2,:) = grad(:,2,:)+reshape(gy,4,1,npts);
grad(:,3,:) = grad(:,3,:)+reshape(gz,4,1,npts);

slp_j = pot(1:3,:);
slp_rho = pot(4,:);
grad_slp_rho = reshape(grad(4,:,:),3,npts);

curl_slp_j = complex(zeros(3,npts));
curl_slp_j(1,:) = reshape(grad(3,2,:)-grad(2,3,:),1,npts);
curl_slp_j(2,:) = reshape(grad(1,3,:)-grad(3,1,:),1,npts);
curl_slp_j(3,:) = reshape(grad(2,1,:)-grad(1,2,:),1,npts);
div_slp_j = reshape( ...
    grad(1,1,:)+grad(2,2,:)+grad(3,3,:),1,npts);

E = 1i*op.zk*slp_j-grad_slp_rho;
nxH = cross(op.geom.n,curl_slp_j,1);
nE = sum(op.geom.n.*E,1);
nxnxE = op.geom.n.*nE-E;

pv = complex(zeros(3,npts));
eq1 = -nxH+op.alpha*nxnxE;
pv(1,:) = sum(op.geom.ru.*eq1,1);
pv(2,:) = sum(op.geom.rv.*eq1,1);
pv(3,:) = -nE+op.alpha*(div_slp_j-1i*op.zk*slp_rho);

y = 0.5*dens_uv+pv;
y = y(:);
end
