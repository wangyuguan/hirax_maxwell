function result = test_ellipsoid_scattering(norder,nproxy)
% TEST_ELLIPSOID_SCATTERING Square sphere-to-sphere Maxwell scattering map.
%
%  Run startup.m first, then result = test_ellipsoid_scattering(16,20).
%  norder controls the physical ellipsoid mesh; nproxy controls the sphere.
%  This small example uses a dense NRCCIE factorization, with no skeletons.

    if nargin < 1, norder = 16; end
    if nargin < 2, nproxy = 20; end
    zk = 2*pi;
    alpha = 1;
    eps_quad = 1e-10;
    axes = [0.5;0.4;0.3];
    radius = 1;
    S = geometries.ellipsoid(axes,[1;1;1],[],norder,11);
    P = maxwell_sphere(radius,nproxy);
    L = maxwell_trace_matrix(P);
    R = nrccie_rhs_matrix(S,alpha);
    fprintf('Ellipsoid: %d nodes, %d NRCCIE unknowns\n',S.npts,3*S.npts);
    fprintf('Sphere: %d nodes, square matrix %d by %d\n', ...
        P.npts,4*P.npts,4*P.npts);

    % Offline: incident sphere trace -> physical RHS -> scattered sphere trace.
    timer = tic;
    corrections = nrccie_quad_corr_block(S,eps_quad,zk,[],128);
    opts = struct('zk',zk,'alpha',alpha,'jump',0.5);
    ids = 1:3*S.npts;
    A = nrccie_system_block(S,[],corrections,opts,ids,ids);
    F = decomposition(A,'lu');
    Qin = -stratton_chu_matrix(P,S,zk);
    B = R*Qin;
    C = L*nrccie_matrix(S,P,zk);
    Scat = C*(F\B);
    build_time = toc(timer);
    fprintf('Built scattering matrix in %.1f s\n',build_time);

    % Three plane waves, including complex polarization, and one exterior dipole.
    % None is used as a training field in the matrix construction.
    labels = {'plane_x','plane_oblique','plane_complex','external_dipole'};
    nsamples = numel(labels);
    incident_proxy = complex(zeros(4*P.npts,nsamples));
    incident_body = complex(zeros(6*S.npts,nsamples));
    for j = 1:nsamples
      incident_proxy(:,j) = L*incident_field(P.r,zk,j);
      incident_body(:,j) = incident_field(S.r,zk,j);
    end

    % Online: the new incident data are sampled only on the enclosing sphere.
    scattered_proxy = Scat*incident_proxy;
    directions = target_directions(80);
    target_radii = [2,4,8];
    T.r = [2*directions,4*directions,8*directions];
    Qout = stratton_chu_matrix(P,T,zk);
    scattered_target = Qout*scattered_proxy;

    % Validation only: solve each actual incident field directly on the ellipsoid.
    rhs = R*incident_body;
    density = F\rhs;
    direct_proxy = C*density;
    direct_target = nrccie_matrix(S,T,zk)*density;
    reconstructed_target = Qout*direct_proxy;
    errors.incoming = relative_errors(Qin*incident_proxy,incident_body);
    errors.trace = relative_errors(scattered_proxy,direct_proxy);
    errors.target = relative_errors(scattered_target,direct_target);
    errors.stratton_chu = relative_errors(reconstructed_target,direct_target);
    errors.bie_residual = relative_errors(A*density,rhs);
    errors.target_e = relative_errors(field_part(scattered_target,1), ...
        field_part(direct_target,1));
    errors.target_h = relative_errors(field_part(scattered_target,2), ...
        field_part(direct_target,2));

    % Analytic outgoing dipole checks Stratton-Chu independently of the BIE.
    source.r = [0.1;-0.1;0.05];
    source.edips = [1;0.3i;-0.2];
    [E,H] = em3d.incoming_sources(zk,source,P,'ed');
    trace = L*reshape([E;H],[],1);
    [E,H] = em3d.incoming_sources(zk,source,T,'ed');
    errors.dipole_identity = relative_errors(Qout*trace,reshape([E;H],[],1));

    fprintf('                       incoming       trace       target       SC only\n');
    for j = 1:nsamples
      fprintf('%-20s  %.3e    %.3e    %.3e    %.3e\n',labels{j}, ...
          errors.incoming(j),errors.trace(j),errors.target(j), ...
          errors.stratton_chu(j));
    end
    fprintf('Analytic outgoing dipole identity: %.3e\n',errors.dipole_identity);
    fprintf('Maximum relative BIE residual: %.3e\n',max(errors.bie_residual));
    assert(all(isfinite(Scat(:))),'Nonfinite scattering matrix.');
    assert(max(errors.bie_residual)<1e-10,'NRCCIE solve failed.');
    if nproxy >= 16
      assert(max(errors.incoming)<2e-4,'Incoming Green reconstruction failed.');
      assert(max(errors.trace)<1e-6,'Scattered sphere trace check failed.');
      assert(errors.dipole_identity<1e-6,'Outgoing Green reconstruction failed.');
    end
    if norder >= 16 && nproxy >= 20
      assert(max(errors.target)<1e-6,'End-to-end target check failed.');
    end

    result = struct('Scat',Scat,'proxy',P,'zk',zk,'alpha',alpha, ...
        'semiaxes',axes,'norder',norder,'nproxy',nproxy, ...
        'physical_nodes',S.npts,'build_time',build_time,'labels',{labels}, ...
        'incident_proxy',incident_proxy,'scattered_proxy',scattered_proxy, ...
        'targets',T,'target_radii',target_radii, ...
        'scattered_target',scattered_target,'direct_target',direct_target, ...
        'errors',errors);
    output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    filename = fullfile(output_dir,sprintf('ellipsoid_p%d_q%d.mat',norder,nproxy));
    save(filename,'-struct','result','-v7.3');
    fprintf('Saved %s\n',filename);
end

function field = incident_field(r,zk,index)
    if index==4
      source.r = [2.5;0.4;-0.3];
      source.edips = [0.2;1;0.3i];
      [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'ed');
    else
      directions = [1,1,-2;0,2,1;0,3,2];
      direction = directions(:,index);
      direction = direction/norm(direction);
      polarization = [0;1;0];
      if index==3, polarization = [1;1i;0]; end
      polarization = polarization-direction*(direction.'*polarization);
      polarization = polarization/norm(polarization);
      phase = exp(1i*zk*(direction.'*r));
      E = polarization*phase;
      H = cross(direction,polarization)*phase;
    end
    field = reshape([E;H],[],1);
end

function directions = target_directions(n)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5));
    s = sqrt(1-z.^2);
    directions = [s.*cos(phi);s.*sin(phi);z];
end

function error = relative_errors(value,reference)
    error = vecnorm(value-reference,2,1)./vecnorm(reference,2,1);
end

function part = field_part(field,index)
    field = reshape(field,6,[],size(field,2));
    part = reshape(field(3*index-2:3*index,:,:),[],size(field,3));
end
