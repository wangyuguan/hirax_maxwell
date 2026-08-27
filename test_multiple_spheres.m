% Test PEC NRCCIE self-convergence on two close spheres.

clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath('../FMM3D/matlab')

% geometry, refinement levels, and tolerances

a = 1;
d = a/30;
na = [4;8;8];
norders = [4,6,8];
iptype = 11;
tol_geom = 100*eps(a);

c = a+d/2;
c0 = [-c,c;0,0;0,0];

wavelength = 2*a;
zk = 2*pi/wavelength;
alpha = 1;
eps = 1e-12;
eps_fmm = 1e-12;
eps_gmres = 1e-11;
maxit = 100;
tol_rres = 5e-11;

abc = a*ones(3,1);
dcheck = norm(c0(:,2)-c0(:,1))-2*a;
assert(dcheck > 0, 'The two spheres overlap or touch.')
assert(abs(dcheck-d) < tol_geom, ...
    'The constructed surface gap differs from the requested gap.')

% Exterior electric point source. NRCCIE requires the incident field to be
% source-free inside every closed component.
src_info = [];
src_info.r = a*[0.30;-0.40;2.80];
src_info.edips = [1+0.30i;-0.40+0.20i;0.25-0.55i];
src_clearance = min(vecnorm(src_info.r-c0,2,1)-a);
assert(src_clearance > 0, ...
    'The incident point source must lie outside both spheres.')

% Fixed exterior targets used to compare the scattered fields.
xyz_out = a*[ ...
    -3.20, 3.20, 0.00,  0.00, 0.00,  0.00; ...
     0.00, 0.00, 3.00, -3.00, 0.00,  0.00; ...
     0.00, 0.00, 0.00,  0.00, 3.00, -3.00];
dist1 = vecnorm(xyz_out-c0(:,1),2,1)-a;
dist2 = vecnorm(xyz_out-c0(:,2),2,1)-a;
target_clearance = min([dist1,dist2]);
assert(target_clearance > 0, ...
    'Every convergence target must lie strictly outside both spheres.')

nlevels = numel(norders);
assert(nlevels >= 3 && all(diff(norders) > 0), ...
    'Use at least three strictly increasing surface orders.')

fields_scattered = complex(zeros(6,size(xyz_out,2),nlevels));
npoints = zeros(1,nlevels);
niters = zeros(1,nlevels);
rresiduals = zeros(1,nlevels);
quadrature_times = zeros(1,nlevels);
solve_times = zeros(1,nlevels);

fprintf('PEC NRCCIE two-sphere self-convergence test\n')
fprintf('  radius: %.8g, gap R/30: %.8g\n',a,dcheck)
fprintf('  wavelength: %.8g = 2R, zk: %.8g\n',wavelength,zk)
fprintf('  source: [%.8g %.8g %.8g], clearance: %.8g\n', ...
    src_info.r,src_clearance)
fprintf('  exterior target clearance: %.8g\n',target_clearance)
fprintf('  subdivisions: [%d %d %d], orders: [%s], iptype: %d\n', ...
    na,num2str(norders),iptype)

for ilevel = 1:nlevels
    norder = norders(ilevel);

    % Extra y/z subdivisions resolve the caps facing across the small gap.
    S1 = geometries.ellipsoid(abc,na,c0(:,1),norder,iptype);
    S2 = geometries.ellipsoid(abc,na,c0(:,2),norder,iptype);
    S = merge([S1,S2]);
    npoints(ilevel) = S.npts;

    n = S.n;
    ru = S.du./vecnorm(S.du,2,1);
    rv = cross(n,ru,1);

    fprintf('\nLevel %d/%d: order %d, patches: %d, points: %d\n', ...
        ilevel,nlevels,norder,S.npatches,S.npts)

    % Add-subtract quadrature corrections.
    rsc = getnear(S);
    itarg = repelem((1:S.npts).',double(diff(rsc.row_ptr)));
    ipatch = double(rsc.col_ind);
    nnear12 = nnz(itarg <= S1.npts & ipatch > S1.npatches);
    nnear21 = nnz(itarg > S1.npts & ipatch <= S1.npatches);
    rfac = rsc.rfac;
    rfac0 = rsc.rfac0;
    clear rsc itarg ipatch

    opts_quad = [];
    opts_quad.rfac = rfac;
    opts_quad.rfac0 = rfac0;

    tquad = tic;
    % Cslp/Cx/Cy/Cz contain accurate-near minus smooth-near.
    Cslp = em3d.slp.get_quad_corr_mat(S,eps,zk,opts_quad);
    [Cx,Cy,Cz] = em3d.sgrad.get_quad_corr_mat(S,eps,zk,opts_quad);
    quadrature_times(ilevel) = toc(tquad);

    i1 = 1:S1.npts;
    i2 = S1.npts+(1:S2.npts);
    ncorr12 = nnz(Cslp(i1,i2))+nnz(Cx(i1,i2))+ ...
        nnz(Cy(i1,i2))+nnz(Cz(i1,i2));
    ncorr21 = nnz(Cslp(i2,i1))+nnz(Cx(i2,i1))+ ...
        nnz(Cy(i2,i1))+nnz(Cz(i2,i1));

    fprintf('  near pairs 1<-2 / 2<-1: %d / %d\n',nnear12,nnear21)
    fprintf('  correction entries 1<-2 / 2<-1: %d / %d\n', ...
        ncorr12,ncorr21)
    fprintf('  correction time: %.2f s\n',quadrature_times(ilevel))
    assert(nnear12 > 0 && nnear21 > 0 && ncorr12 > 0 && ncorr21 > 0, ...
        'Close-surface corrections are missing at order %d.',norder)

    % Incident field and NRCCIE right-hand side.
    [Einc,Hinc] = em3d.incoming_sources( ...
        zk,src_info,S,'electric dipole');
    nEinc = sum(n.*Einc,1);
    nxHinc = cross(n,Hinc,1);
    nxnxEinc = n.*nEinc-Einc;
    rhs_t = nxHinc-alpha*nxnxEinc;

    rhs = complex(zeros(3,S.npts));
    rhs(1,:) = sum(ru.*rhs_t,1);
    rhs(2,:) = sum(rv.*rhs_t,1);
    rhs(3,:) = nEinc;
    rhs = rhs(:);

    % Matrix-free NRCCIE solve.
    op = [];
    op.npts = S.npts;
    op.r = S.r;
    op.wts = S.wts(:).';
    op.n = n;
    op.ru = ru;
    op.rv = rv;
    op.zk = zk;
    op.alpha = alpha;
    op.eps_fmm = eps_fmm;
    op.Cslp = Cslp;
    op.Cx = Cx;
    op.Cy = Cy;
    op.Cz = Cz;

    matvec = @(x) nrccie_matvec(x,op);
    tsolve = tic;
    [soln,flag,~,~,errs] = gmres(matvec,rhs,[],eps_gmres,maxit);
    solve_times(ilevel) = toc(tsolve);

    niters(ilevel) = length(errs)-1;
    rresiduals(ilevel) = norm(matvec(soln)-rhs)/norm(rhs);
    dens_uv = reshape(soln,3,S.npts);
    zjvec = ru.*dens_uv(1,:)+rv.*dens_uv(2,:);
    rho = dens_uv(3,:);
    err_n = max(abs(sum(n.*zjvec,1)));

    fprintf('  GMRES iterations: %d, flag: %d\n',niters(ilevel),flag)
    fprintf('  relative residual: %.3e\n',rresiduals(ilevel))
    fprintf('  maximum normal current: %.3e\n',err_n)
    fprintf('  solve time: %.2f s\n',solve_times(ilevel))

    assert(flag == 0, 'GMRES did not converge at order %d.',norder)
    assert(rresiduals(ilevel) < tol_rres, ...
        'Order-%d residual %.3e exceeds %.3e.', ...
        norder,rresiduals(ilevel),tol_rres)

    [E,H] = eval_nrccie(zjvec,rho,xyz_out,op);
    fields_scattered(:,:,ilevel) = [E;H];

    clear Cslp Cx Cy Cz E H Einc Hinc S S1 S2 matvec op soln
end

% Compare every non-reference level with the finest-order scattered field.
reference_field = fields_scattered(:,:,end);
reference_norm = norm(reference_field,'fro');
assert(reference_norm > 0, 'The reference scattered field is zero.')

self_errors = zeros(1,nlevels-1);
step_changes = zeros(1,nlevels-1);
for ilevel = 1:nlevels-1
    self_errors(ilevel) = norm( ...
        fields_scattered(:,:,ilevel)-reference_field,'fro')/reference_norm;
    step_changes(ilevel) = norm( ...
        fields_scattered(:,:,ilevel+1)-fields_scattered(:,:,ilevel), ...
        'fro')/norm(fields_scattered(:,:,ilevel+1),'fro');
end

fprintf('\nSelf-convergence against order %d\n',norders(end))
for ilevel = 1:nlevels-1
    fprintf('  order %d -> %d: relative field change %.3e\n', ...
        norders(ilevel),norders(ilevel+1),step_changes(ilevel))
    fprintf('  order %d vs %d: relative field error %.3e\n', ...
        norders(ilevel),norders(end),self_errors(ilevel))
end
fprintf('  error reduction from order %d to %d: %.3g\n', ...
    norders(1),norders(end-1),self_errors(1)/self_errors(end))

assert(all(diff(self_errors) < 0), ...
    'Scattered fields do not converge monotonically to the finest level.')
fprintf('  PASS\n')


function y = nrccie_matvec(x,op)

npts = op.npts;
dens_uv = reshape(x,3,npts);
zjvec = op.ru.*dens_uv(1,:)+op.rv.*dens_uv(2,:);
rho = dens_uv(3,:);
dens = [zjvec;rho];

srcinfo = [];
srcinfo.sources = op.r;
srcinfo.nd = 4;
srcinfo.charges = dens.*op.wts;

U = hfmm3d(op.eps_fmm,op.zk,srcinfo,2);
pot = reshape(U.pot,4,npts);
grad = reshape(U.grad,4,3,npts);

pot = pot+(op.Cslp*dens.').';
gx = (op.Cx*dens.').';
gy = (op.Cy*dens.').';
gz = (op.Cz*dens.').';
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
nxH = cross(op.n,curl_slp_j,1);
nE = sum(op.n.*E,1);
nxnxE = op.n.*nE-E;

pv = complex(zeros(3,npts));
eq1 = -nxH+op.alpha*nxnxE;
pv(1,:) = sum(op.ru.*eq1,1);
pv(2,:) = sum(op.rv.*eq1,1);
pv(3,:) = -nE+op.alpha*(div_slp_j-1i*op.zk*slp_rho);

y = 0.5*dens_uv+pv;
y = y(:);
end


function [E,H] = eval_nrccie(zjvec,rho,targs,op)

ntarg = size(targs,2);
dens = [zjvec;rho];

srcinfo = [];
srcinfo.sources = op.r;
srcinfo.nd = 4;
srcinfo.charges = dens.*op.wts;

U = hfmm3d(op.eps_fmm,op.zk,srcinfo,0,targs,2);
pot = reshape(U.pottarg,4,ntarg);
grad = reshape(U.gradtarg,4,3,ntarg);

E = 1i*op.zk*pot(1:3,:)-reshape(grad(4,:,:),3,ntarg);
H = complex(zeros(3,ntarg));
H(1,:) = reshape(grad(3,2,:)-grad(2,3,:),1,ntarg);
H(2,:) = reshape(grad(1,3,:)-grad(3,1,:),1,ntarg);
H(3,:) = reshape(grad(2,1,:)-grad(1,2,:),1,ntarg);
end
