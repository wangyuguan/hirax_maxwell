function result = test_ellipsoid_bc_skeleton(tol,norder,nproxy)
% TEST_ELLIPSOID_BC_SKELETON Compress the sphere interface before BIE solves.
%
%  Input arguments:
%    * tol: relative joint block ID tolerance (1e-6)
%    * norder: physical ellipsoid discretization order (16)
%    * nproxy: Gauss-Legendre order on the enclosing sphere (20)
%
%  Output arguments:
%    * result: square compressed Scat, shared nodes, interface maps, and errors
%
%  Run startup.m first. No full scattering matrix is formed or loaded.
%  The input is input_projection times the full incident trace, not just
%  incident values at the selected nodes. The output contains scattered
%  trace values at those nodes; output_interpolation reconstructs the sphere.

    if nargin < 1 || isempty(tol), tol = 1e-6; end
    if nargin < 2 || isempty(norder), norder = 16; end
    if nargin < 3 || isempty(nproxy), nproxy = 20; end
    zk = 2*pi;
    alpha = 1;
    eps_quad = 1e-10;
    semiaxes = [0.5;0.4;0.3];
    timer = tic;
    S = geometries.ellipsoid(semiaxes,[1;1;1],[],norder,11);
    P = maxwell_sphere(1,nproxy);
    L = maxwell_trace_matrix(P);
    R = nrccie_rhs_matrix(S,alpha);
    fprintf('Ellipsoid: %d nodes, %d NRCCIE unknowns\n',S.npts,3*S.npts);
    fprintf('Sphere: %d nodes, %d trace coordinates\n',P.npts,4*P.npts);

    % B maps incident traces to the physical RHS; C evaluates physical density.
    stage = tic;
    Qin = -stratton_chu_matrix(P,S,zk);
    B = R*Qin;
    clear Qin
    C = L*nrccie_matrix(S,P,zk);
    timings.kernels = toc(stage);
    normB = norm(B,'fro');
    normC = norm(C,'fro');
    fprintf('Built B and C in %.1f s; selecting shared point blocks\n', ...
        timings.kernels);

    % Row ID of [B*,C] gives B approximately B(:,I)*U* and C approximately U*C(I,:).
    stage = tic;
    [nodes,U,id_info] = maxwell_shared_row_id([B'/normB,C/normC],tol);
    dofs = id_info.dofs;
    input_projection = U';
    output_interpolation = U;
    [errors.operator_b,errors.operator_c] = operator_errors(B,C,U,dofs);
    timings.id = toc(stage);
    fprintf('Selected %d/%d nodes; compressed Scat is %d by %d\n', ...
        numel(nodes),P.npts,numel(dofs),numel(dofs));
    fprintf('Operator residuals: B %.3e, C %.3e; ID time %.1f s\n', ...
        errors.operator_b,errors.operator_c,timings.id);

    % Only selected columns of B require BIE solves to construct Scat.
    stage = tic;
    corrections = nrccie_quad_corr_block(S,eps_quad,zk,[],128);
    opts = struct('zk',zk,'alpha',alpha,'jump',0.5);
    ids = 1:3*S.npts;
    A = nrccie_system_block(S,[],corrections,opts,ids,ids);
    F = decomposition(A,'lu');
    clear corrections
    timings.factorization = toc(stage);
    stage = tic;
    Scat = complex(zeros(numel(dofs),numel(dofs)));
    residual2 = 0;
    rhs2 = 0;
    for first = 1:128:numel(dofs)
      cols = first:min(first+127,numel(dofs));
      rhs = B(:,dofs(cols));
      density = F\rhs;
      Scat(:,cols) = C(dofs,:)*density;
      residual2 = residual2+norm(A*density-rhs,'fro')^2;
      rhs2 = rhs2+norm(rhs,'fro')^2;
    end
    errors.bie_residual_small = sqrt(residual2/rhs2);
    timings.selected_solves = toc(stage);
    build_time = toc(timer);
    fprintf('Factorization %.1f s; %d selected RHS solves %.1f s\n', ...
        timings.factorization,numel(dofs),timings.selected_solves);

    % These fields were not used for selecting nodes. Both surfaces are sampled
    % analytically, so the reference does not rely on an incoming Green map.
    stage = tic;
    [incident_proxy,labels] = incident_fields(P.r,zk);
    incident_proxy = L*incident_proxy;
    incident_body = incident_fields(S.r,zk);
    rhs = R*incident_body;
    rhs_green = B*incident_proxy;
    incident_coefficients = input_projection*incident_proxy;
    rhs_compressed = B(:,dofs)*incident_coefficients;
    density = F\rhs;
    density_green = F\rhs_green;
    direct_proxy = C*density;
    reference_proxy = C*density_green;
    scattered_skeleton = Scat*incident_coefficients;
    scattered_proxy = output_interpolation*scattered_skeleton;

    target_radii = [2,4,8];
    d = directions(80,0.31);
    T.r = [2*d,4*d,8*d];
    Q = stratton_chu_matrix(P,T,zk);
    target_matrix = Q*output_interpolation;
    scattered_target = target_matrix*scattered_skeleton;
    reference_target = Q*reference_proxy;
    direct_target = nrccie_matrix(S,T,zk)*density;
    errors.incoming_green_rhs = relative_errors(rhs_green,rhs);
    errors.incoming_rhs_compression = relative_errors(rhs_compressed,rhs_green);
    errors.outgoing_trace = relative_errors(U*direct_proxy(dofs,:),direct_proxy);
    errors.skeleton = relative_errors(scattered_skeleton,reference_proxy(dofs,:));
    errors.proxy_compression = relative_errors(scattered_proxy,reference_proxy);
    errors.proxy_direct_bie = relative_errors(scattered_proxy,direct_proxy);
    errors.target_compression = relative_errors(scattered_target,reference_target);
    errors.target_direct_bie = relative_errors(scattered_target,direct_target);
    errors.baseline_target_direct_bie = relative_errors(reference_target,direct_target);
    errors.bie_residual = relative_errors(A*density,rhs);
    errors.bie_residual_green = relative_errors(A*density_green,rhs_green);
    timings.validation = toc(stage);
    fprintf('Validation: %d independent fields, %.1f s\n',numel(labels),timings.validation);
    fprintf('Incoming Green RHS: %.3e; incoming compression: %.3e\n', ...
        max(errors.incoming_green_rhs),max(errors.incoming_rhs_compression));
    fprintf('Outgoing trace interpolation: %.3e; combined trace: %.3e\n', ...
        max(errors.outgoing_trace),max(errors.proxy_compression));
    fprintf('Added target error %.3e; total target error %.3e; baseline %.3e\n', ...
        max(errors.target_compression),max(errors.target_direct_bie), ...
        max(errors.baseline_target_direct_bie));

    skeleton = struct('r',P.r(:,nodes),'n',P.n(:,nodes), ...
        'dru',P.dru(:,nodes),'drv',P.drv(:,nodes),'npts',numel(nodes));
    normalization = struct('b_frobenius',normB,'c_frobenius',normC);
    result = struct('Scat',Scat,'nodes',nodes,'dofs',dofs, ...
        'input_projection',input_projection,'output_interpolation',output_interpolation, ...
        'skeleton',skeleton,'proxy',P,'zk',zk,'alpha',alpha,'tol',tol, ...
        'semiaxes',semiaxes,'physical_nodes',S.npts,'norder',norder,'nproxy',nproxy, ...
        'id_info',id_info,'normalization',normalization,'build_time',build_time, ...
        'timings',timings,'labels',{labels},'errors',errors, ...
        'incident_proxy',incident_proxy,'incident_coefficients',incident_coefficients, ...
        'scattered_skeleton',scattered_skeleton,'scattered_proxy',scattered_proxy, ...
        'targets',T,'target_radii',target_radii,'target_matrix',target_matrix, ...
        'scattered_target',scattered_target,'reference_target',reference_target, ...
        'direct_target',direct_target);
    output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    tag = strrep(sprintf('%.0e',tol),'-','m');
    filename = fullfile(output_dir,['ellipsoid_bc_skeleton_',tag,'.mat']);
    save(filename,'-struct','result','-v7.3');
    fprintf('Saved %s\n',filename);

    assert(isequal(size(Scat),[numel(dofs),numel(dofs)]));
    assert(norm(U(dofs,:)-eye(numel(dofs)),'fro')<1e-12);
    assert(all(isfinite(Scat(:))),'Nonfinite compressed scattering matrix.');
    assert(id_info.relative_residual<tol,'Joint operator ID failed.');
    assert(max([errors.operator_b,errors.operator_c])<2*tol,'Operator ID residual failed.');
    assert(numel(nodes)<P.npts,'No reduction in sphere node count.');
    assert(errors.bie_residual_small<1e-10,'Selected-column BIE solves failed.');
    assert(max(errors.bie_residual)<1e-10,'Independent BIE solves failed.');
    assert(max(errors.bie_residual_green)<1e-10,'Green-input BIE solves failed.');
    assert(max(errors.incoming_rhs_compression)<100*tol,'Incoming RHS compression failed.');
    assert(max(errors.outgoing_trace)<100*tol,'Outgoing trace interpolation failed.');
    assert(max(errors.proxy_compression)<100*tol,'Compressed sphere trace check failed.');
    assert(max(errors.target_compression)<100*tol,'Target compression check failed.');
    if norder>=16 && nproxy>=20
      assert(max(errors.target_direct_bie)<2e-6,'End-to-end target check failed.');
    end
end

function [errB,errC] = operator_errors(B,C,U,dofs)
    residualB = 0;
    residualC = 0;
    for first = 1:128:size(U,1)
      rows = first:min(first+127,size(U,1));
      residualB = residualB+norm(B(:,rows)-B(:,dofs)*U(rows,:)','fro')^2;
      residualC = residualC+norm(C(rows,:)-U(rows,:)*C(dofs,:),'fro')^2;
    end
    errB = sqrt(residualB)/norm(B,'fro');
    errC = sqrt(residualC)/norm(C,'fro');
end

function [fields,labels] = incident_fields(r,zk)
    d = directions(16,0.73);
    fields = complex(zeros(6*size(r,2),48));
    labels = cell(1,48);
    for j = 1:16
      [~,axis] = min(abs(d(:,j)));
      unit = zeros(3,1); unit(axis) = 1;
      u = cross(d(:,j),unit); u = u/norm(u);
      v = cross(d(:,j),u);
      phase = exp(1i*zk*(d(:,j).'*r));
      fields(:,2*j-1) = reshape([u*phase;v*phase],[],1);
      fields(:,2*j) = reshape([v*phase;-u*phase],[],1);
      labels{2*j-1} = sprintf('plane_%02d_u',j);
      labels{2*j} = sprintf('plane_%02d_v',j);
    end
    for j = 1:4
      fields(:,32+j) = (fields(:,2*j-1)+1i*fields(:,2*j))/sqrt(2);
      labels{32+j} = sprintf('complex_plane_%02d',j);
    end
    positions = directions(2,0.61).*[2.5,3];
    for j = 1:2
      source.r = positions(:,j);
      for c = 1:3
        dipole = zeros(3,1); dipole(c) = 1;
        source.edips = dipole;
        source.hdips = dipole;
        [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'ed');
        fields(:,36+6*j-6+c) = reshape([E;H],[],1);
        labels{36+6*j-6+c} = sprintf('electric_dipole_%d_%d',j,c);
        [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'hd');
        fields(:,36+6*j-3+c) = reshape([E;H],[],1);
        labels{36+6*j-3+c} = sprintf('magnetic_dipole_%d_%d',j,c);
      end
    end
end

function d = directions(n,angle)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5))+angle;
    s = sqrt(1-z.^2);
    d = [s.*cos(phi);s.*sin(phi);z];
    rotation = [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
    d = rotation*d;
end

function error = relative_errors(value,reference)
    error = vecnorm(value-reference,2,1)./vecnorm(reference,2,1);
end
