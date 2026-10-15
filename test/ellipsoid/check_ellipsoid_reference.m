function diagnostics = check_ellipsoid_reference(filename,norder)
%
%  check_ellipsoid_reference
%    Compare saved ellipsoid fields with the independent fmm3dbie solver.
%
%  Syntax
%    diagnostics = check_ellipsoid_reference(filename,norder)
%
%  Input arguments:
%    * filename: saved ellipsoid_p<N>_q<M>.mat from the scattering test
%        (test/results/ellipsoid_p16_q20.mat)
%    * norder: reference ellipsoid order (max(16,saved order))
%
%  Output arguments:
%    * diagnostics: field errors, reference fields, residual and timings
%
%  Run startup.m first. The first saved incident field has direction +x,
%  electric polarization +y and magnetic polarization +z. This check uses
%  em3d.pec.solver/eval, including their oversampled quadrature, rather than
%  the dense system and field matrices used by the scattering test. Agreement
%  at one order checks the implementations; mesh refinement is still needed
%  to establish accuracy for the continuous scattering problem.

    output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if nargin < 1 || isempty(filename)
      filename = fullfile(output_dir,'ellipsoid_p16_q20.mat');
    end
    saved = load(filename,'zk','alpha','semiaxes','norder','nproxy', ...
        'labels','targets','target_radii','direct_target','scattered_target');
    assert(strcmp(saved.labels{1},'plane_x'), ...
        'The first saved incident field must be the +x plane wave.');
    if nargin < 2 || isempty(norder)
      norder = max(16,saved.norder);
    end
    eps_quad = 1e-10;
    opts = struct('eps_gmres',1e-11,'maxit',200,'rep','nrccie');
    S = geometries.ellipsoid(saved.semiaxes,[1;1;1],[],norder,11);
    phase = exp(1i*saved.zk*S.r(1,:));
    Einc = [0;1;0]*phase;
    Hinc = [0;0;1]*phase;
    fprintf('Independent reference: order %d, %d physical nodes\n', ...
        norder,S.npts);

    timer = tic;
    [density,history,residual] = em3d.pec.solver( ...
        S,Einc,Hinc,eps_quad,saved.zk,saved.alpha,opts);
    timings.solve = toc(timer);
    timer = tic;
    [E,H] = em3d.pec.eval(S,density,saved.targets, ...
        eps_quad,saved.zk,saved.alpha);
    timings.evaluate = toc(timer);
    reference = reshape([E;H],[],1);
    direct = reshape(saved.direct_target(:,1),6,[]);
    scattered = reshape(saved.scattered_target(:,1),6,[]);
    errors.direct = relative_error(direct,[E;H]);
    errors.scattering = relative_error(scattered,[E;H]);
    errors.direct_e = relative_error(direct(1:3,:),E);
    errors.direct_h = relative_error(direct(4:6,:),H);
    errors.scattering_e = relative_error(scattered(1:3,:),E);
    errors.scattering_h = relative_error(scattered(4:6,:),H);
    radii = vecnorm(saved.targets.r,2,1);
    errors.by_radius = zeros(numel(saved.target_radii),4);
    for j = 1:numel(saved.target_radii)
      radius = saved.target_radii(j);
      ids = abs(radii-radius)<1e-12*max(1,radius);
      assert(any(ids),'No targets found at a saved target radius.');
      errors.by_radius(j,:) = [ ...
          relative_error(direct(1:3,ids),E(:,ids)), ...
          relative_error(direct(4:6,ids),H(:,ids)), ...
          relative_error(scattered(1:3,ids),E(:,ids)), ...
          relative_error(scattered(4:6,ids),H(:,ids))];
    end
    fprintf('Reference GMRES: %d iterations, residual %.3e\n', ...
        numel(history),residual);
    fprintf('Compared with reference:       E             H\n');
    fprintf('Direct physical evaluation  %.3e     %.3e\n', ...
        errors.direct_e,errors.direct_h);
    fprintf('Square map + Stratton-Chu   %.3e     %.3e\n', ...
        errors.scattering_e,errors.scattering_h);
    fprintf('Radius        direct E      direct H      square E      square H\n');
    for j = 1:numel(saved.target_radii)
      fprintf('%6.2f        %.3e     %.3e     %.3e     %.3e\n', ...
          saved.target_radii(j),errors.by_radius(j,:));
    end

    diagnostics = struct('input_file',filename,'norder',norder, ...
        'physical_nodes',S.npts,'saved_norder',saved.norder, ...
        'saved_nproxy',saved.nproxy,'eps_quad',eps_quad, ...
        'gmres_tolerance',opts.eps_gmres,'gmres_history',history, ...
        'gmres_residual',residual,'reference_target',reference, ...
        'target_radii',saved.target_radii, ...
        'radius_error_columns',{{'direct_e','direct_h','scattering_e','scattering_h'}}, ...
        'errors',errors,'timings',timings);
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    [~,name] = fileparts(filename);
    output_file = fullfile(output_dir,sprintf('%s_reference_p%d.mat',name,norder));
    save(output_file,'-struct','diagnostics');
    fprintf('Saved %s\n',output_file);
    assert(all(isfinite(reference)) && isfinite(residual), ...
        'The independent reference produced nonfinite values.');
    assert(residual<10*opts.eps_gmres,'The independent GMRES reference did not converge.');
end

function value = relative_error(actual,reference)
    value = norm(actual-reference,'fro')/max(norm(reference,'fro'),realmin);
end
