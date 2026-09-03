% Compare original and rotationally reused four-object NRCCIE systems.

clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath('../FMM3D/matlab')
addpath('src')

rng(1)

a = 1;
gap = a/30;
center_offset = a+gap/2;
norder = 2;
sphere_subdivisions = 1;
iptype = 1;
zk = 2*pi/(10*a);
alpha = 1;
eps_quad = 1e-8;
eps_fmm = 1e-12;
eps_gmres = 1e-10;
maxit = 200;

S0 = geometries.sphere( ...
    a,sphere_subdivisions,[0;0;0],norder,iptype);

theta = [pi/4,-pi/4,3*pi/4,-3*pi/4];
shifts = center_offset*[ ...
    -1, 1,-1, 1; ...
     1, 1,-1,-1; ...
     0, 0, 0, 0];
surfaces = cell(1,4);
for isurface = 1:4
    ct = cos(theta(isurface));
    st = sin(theta(isurface));
    rotation = [ct,-st,0;st,ct,0;0,0,1];
    surfaces{isurface} = affine_transf( ...
        S0,rotation,shifts(:,isurface));
end
S1 = surfaces{1};
S2 = surfaces{2};
S3 = surfaces{3};
S4 = surfaces{4};
S = merge([S1,S2,S3,S4]);

fprintf('Four-sphere correction-symmetry A/B test\n')
fprintf('  nodes per sphere / total: %d / %d\n',S0.npts,S.npts)

% Original full-surface construction.
tfull = tic;
Cslp_full = em3d.slp.get_quad_corr_mat(S,eps_quad,zk);
[Cx_full,Cy_full,Cz_full] = ...
    em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk);
tfull = toc(tfull);

% One target block row followed by rotational reuse.
tsym = tic;
[Cslp_sym,Cx_sym,Cy_sym,Cz_sym] = ...
    fourfold_quad_corr_mats(S1,S2,S3,S4,eps_quad,zk);
tsym = toc(tsym);

C_full = {Cslp_full,Cx_full,Cy_full,Cz_full};
C_sym = {Cslp_sym,Cx_sym,Cy_sym,Cz_sym};
correction_errors = zeros(1,4);
for ikernel = 1:4
    correction_errors(ikernel) = norm( ...
        C_sym{ikernel}-C_full{ikernel},'fro')/ ...
        norm(C_full{ikernel},'fro');
end

fprintf('  full / reused construction: %.3f / %.3f s\n',tfull,tsym)
fprintf('  relative correction errors S / dx / dy / dz:\n')
fprintf('    %.3e / %.3e / %.3e / %.3e\n',correction_errors)

surface_normal = S.n;
surface_ru = S.du./vecnorm(S.du,2,1);
surface_rv = cross(surface_normal,surface_ru,1);

operator = struct();
operator.npts = S.npts;
operator.r = S.r;
operator.wts = S.wts(:).';
operator.n = surface_normal;
operator.ru = surface_ru;
operator.rv = surface_rv;
operator.zk = zk;
operator.alpha = alpha;
operator.eps_fmm = eps_fmm;

operator_full = operator;
operator_full.C = C_full;
operator_sym = operator;
operator_sym.C = C_sym;

matvec_full = @(x) apply_nrccie_test(x,operator_full);
matvec_sym = @(x) apply_nrccie_test(x,operator_sym);

% This density is arbitrary and has no rotational relationship between
% objects, so the operator comparison does not use data symmetry.
x_test = randn(3*S.npts,1)+1i*randn(3*S.npts,1);
y_full = matvec_full(x_test);
y_sym = matvec_sym(x_test);
matvec_error = norm(y_sym-y_full)/norm(y_full);
fprintf('  arbitrary-density matvec error: %.3e\n',matvec_error)

% Use an off-axis exterior source so the right-hand side is not invariant
% under the fourfold geometry rotation.
source_info = struct();
source_info.r = a*[0.30;-0.40;2.80];
source_info.edips = [1+0.30i;-0.40+0.20i;0.25-0.55i];
[einc,hinc] = em3d.incoming_sources( ...
    zk,source_info,S,'electric dipole');

normal_einc = sum(surface_normal.*einc,1);
nxhinc = cross(surface_normal,hinc,1);
nxnxeinc = surface_normal.*normal_einc-einc;
tangent_rhs = nxhinc-alpha*nxnxeinc;
rhs_components = complex(zeros(3,S.npts));
rhs_components(1,:) = sum(surface_ru.*tangent_rhs,1);
rhs_components(2,:) = sum(surface_rv.*tangent_rhs,1);
rhs_components(3,:) = normal_einc;
rhs = rhs_components(:);

[solution_full,flag_full] = gmres( ...
    matvec_full,rhs,[],eps_gmres,maxit);
[solution_sym,flag_sym] = gmres( ...
    matvec_sym,rhs,[],eps_gmres,maxit);

solution_error = norm(solution_sym-solution_full)/norm(solution_full);
full_residual = norm(matvec_full(solution_full)-rhs)/norm(rhs);
sym_residual_in_full_system = ...
    norm(matvec_full(solution_sym)-rhs)/norm(rhs);

fprintf('  GMRES flags original / reused: %d / %d\n', ...
    flag_full,flag_sym)
fprintf('  solution difference: %.3e\n',solution_error)
fprintf('  original-system residual, original solution: %.3e\n', ...
    full_residual)
fprintf('  original-system residual, reused solution: %.3e\n', ...
    sym_residual_in_full_system)

assert(max(correction_errors) < 1e-10, ...
    'Rotationally reused correction matrices differ from the originals.')
assert(matvec_error < 1e-10, ...
    'The complete NRCCIE operators differ.')
assert(flag_full == 0 && flag_sym == 0, ...
    'One of the two GMRES solves did not converge.')
assert(solution_error < 1e-8, ...
    'The original and reused systems produce different solutions.')
assert(sym_residual_in_full_system < 1e-9, ...
    'The reused solution does not satisfy the original system.')
fprintf('  PASS\n')


function y = apply_nrccie_test(x,operator)

npts = operator.npts;
density_components = reshape(x,3,npts);
surface_current = operator.ru.*density_components(1,:)+ ...
    operator.rv.*density_components(2,:);
surface_charge = density_components(3,:);
density = [surface_current;surface_charge];

source = struct();
source.sources = operator.r;
source.nd = 4;
source.charges = density.*operator.wts;

fmm_output = hfmm3d(operator.eps_fmm,operator.zk,source,2);
potential = reshape(fmm_output.pot,4,npts);
gradient = reshape(fmm_output.grad,4,3,npts);

potential = potential+(operator.C{1}*density.').';
for idirection = 1:3
    correction = (operator.C{idirection+1}*density.').';
    gradient(:,idirection,:) = gradient(:,idirection,:)+ ...
        reshape(correction,4,1,npts);
end

slp_current = potential(1:3,:);
slp_charge = potential(4,:);
gradient_slp_charge = reshape(gradient(4,:,:),3,npts);

curl_slp_current = complex(zeros(3,npts));
curl_slp_current(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,npts);
curl_slp_current(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,npts);
curl_slp_current(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,npts);
divergence_slp_current = reshape( ...
    gradient(1,1,:)+gradient(2,2,:)+gradient(3,3,:),1,npts);

electric_field = 1i*operator.zk*slp_current-gradient_slp_charge;
nxh = cross(operator.n,curl_slp_current,1);
normal_electric_field = sum(operator.n.*electric_field,1);
nxnxe = operator.n.*normal_electric_field-electric_field;

principal_value = complex(zeros(3,npts));
tangent_equation = -nxh+operator.alpha*nxnxe;
principal_value(1,:) = sum(operator.ru.*tangent_equation,1);
principal_value(2,:) = sum(operator.rv.*tangent_equation,1);
principal_value(3,:) = -normal_electric_field+operator.alpha*( ...
    divergence_slp_current-1i*operator.zk*slp_charge);

y = 0.5*density_components+principal_value;
y = y(:);
end
