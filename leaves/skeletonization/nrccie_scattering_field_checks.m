function [checks,fields] = nrccie_scattering_field_checks(S,outgoing_source, ...
    zk,center,directions,radii,reference_density,direct_density, ...
    equivalent_density)
%NRCCIE_SCATTERING_FIELD_CHECKS Compare scattering and full-density E/H fields.
%
% Input arguments:
%   S, outgoing_source - full surface and outgoing skeleton nodes.
%   zk - wave number.
%   center, directions, radii - targets are center + radius*directions,
%                              with unit columns in directions.
%   reference_density - plain GMRES density; [] skips reference comparisons.
%   direct_density - full-surface direct RS-S density.
%   equivalent_density - scattering-matrix output on outgoing_source.
%
% Output arguments:
%   checks - separate E/H relative errors, with size nradii-by-nrhs.
%            A_vs_B denotes norm(B-A)/max(norm(A),realmin).
%   fields - reference, direct and scattering fields, each of size
%            6-by-ndirections-by-nrhs-by-nradii, with components
%            [Ex;Ey;Ez;Hx;Hy;Hz]. Disabled reference arrays are empty.
%
% Fields are evaluated directly at finite radii, without radial rescaling.

    validateattributes(zk,{'numeric'},{'scalar','finite','nonzero'})
    validateattributes(center,{'numeric'},{'vector','numel',3,'real','finite'})
    validateattributes(directions,{'numeric'}, ...
      {'2d','nrows',3,'nonempty','real','finite'})
    validateattributes(radii,{'numeric'}, ...
      {'vector','nonempty','real','finite','positive'})
    center = reshape(center,3,1);
    radii = reshape(radii,1,[]);
    direction_norms = vecnorm(directions,2,1);
    assert(all(abs(direction_norms-1)<=1e-12), ...
      'nrccie_scattering_field_checks:NonunitDirections', ...
      'Each direction must have unit Euclidean norm.')
    assert(size(direct_density,1)==3*S.npts && ...
      size(equivalent_density,1)==3*outgoing_source.npts, ...
      'Each density must contain three interleaved rows per source node.')
    reference_enabled = ~isempty(reference_density);
    if reference_enabled
      assert(size(reference_density,1)==3*S.npts, ...
        'The reference density must contain three rows per full-surface node.')
    end
    nwave = size(direct_density,2);
    assert(nwave>0 && size(equivalent_density,2)==nwave && ...
      (~reference_enabled || size(reference_density,2)==nwave), ...
      'nrccie_scattering_field_checks:InconsistentRHS', ...
      'All enabled densities must contain the same positive number of RHS columns.')
    ndir = size(directions,2);
    nradii = numel(radii);
    reference_errors = [];
    if reference_enabled
      reference_errors = zeros(nradii,nwave);
    end

    checks = struct('center',center,'directions',directions,'radii',radii, ...
      'zk',zk,'wavelength',NaN,'radii_in_wavelengths',NaN(size(radii)), ...
      'reference_enabled',reference_enabled, ...
      'reference_vs_direct_e',reference_errors, ...
      'reference_vs_direct_h',reference_errors, ...
      'reference_vs_scattering_e',reference_errors, ...
      'reference_vs_scattering_h',reference_errors, ...
      'direct_vs_scattering_e',zeros(nradii,nwave), ...
      'direct_vs_scattering_h',zeros(nradii,nwave), ...
      'radius_seconds',zeros(nradii,1),'seconds',0, ...
      'finite_radius',true,'asymptotic_extrapolation',false);
    if isreal(zk) && zk>0
      checks.wavelength = 2*pi/zk;
      checks.radii_in_wavelengths = radii/checks.wavelength;
    end
    fields = struct();
    fields.reference = [];
    if reference_enabled
      fields.reference = complex(zeros(6,ndir,nwave,nradii));
    end
    fields.direct = complex(zeros(6,ndir,nwave,nradii));
    fields.scattering = complex(zeros(6,ndir,nwave,nradii));
    fields.targets = zeros(3,ndir,nradii);
    total_timer = tic;
    for ir = 1:nradii
      radius_timer = tic;
      targets = struct('r',center+radii(ir)*directions);
      fields.targets(:,:,ir) = targets.r;
      Kfull = nrccie_matrix(S,targets,zk);
      Kout = nrccie_matrix(outgoing_source,targets,zk);
      if reference_enabled
        fields.reference(:,:,:,ir) = reshape(Kfull*reference_density,6,ndir,nwave);
      end
      fields.direct(:,:,:,ir) = reshape(Kfull*direct_density,6,ndir,nwave);
      fields.scattering(:,:,:,ir) = reshape(Kout*equivalent_density,6,ndir,nwave);
      clear Kfull Kout
      for wave = 1:nwave
        direct_field = fields.direct(:,:,wave,ir);
        scattering = fields.scattering(:,:,wave,ir);
        if reference_enabled
          reference = fields.reference(:,:,wave,ir);
          checks.reference_vs_direct_e(ir,wave) = ...
            relative_error(direct_field(1:3,:),reference(1:3,:));
          checks.reference_vs_direct_h(ir,wave) = ...
            relative_error(direct_field(4:6,:),reference(4:6,:));
          checks.reference_vs_scattering_e(ir,wave) = ...
            relative_error(scattering(1:3,:),reference(1:3,:));
          checks.reference_vs_scattering_h(ir,wave) = ...
            relative_error(scattering(4:6,:),reference(4:6,:));
        end
        checks.direct_vs_scattering_e(ir,wave) = ...
          relative_error(scattering(1:3,:),direct_field(1:3,:));
        checks.direct_vs_scattering_h(ir,wave) = ...
          relative_error(scattering(4:6,:),direct_field(4:6,:));
      end
      checks.radius_seconds(ir) = toc(radius_timer);
      if reference_enabled
        fprintf(['Finite-radius E/H at r=%.8g: direct/reference %.3e/%.3e; ' ...
          'Scat/reference %.3e/%.3e; Scat/direct %.3e/%.3e (%.2f s).\n'], ...
          radii(ir),max(checks.reference_vs_direct_e(ir,:)), ...
          max(checks.reference_vs_direct_h(ir,:)), ...
          max(checks.reference_vs_scattering_e(ir,:)), ...
          max(checks.reference_vs_scattering_h(ir,:)), ...
          max(checks.direct_vs_scattering_e(ir,:)), ...
          max(checks.direct_vs_scattering_h(ir,:)),checks.radius_seconds(ir))
      else
        fprintf(['Finite-radius E/H at r=%.8g: Scat/direct %.3e/%.3e ' ...
          '(reference disabled; %.2f s).\n'],radii(ir), ...
          max(checks.direct_vs_scattering_e(ir,:)), ...
          max(checks.direct_vs_scattering_h(ir,:)),checks.radius_seconds(ir))
      end
    end
    checks.seconds = toc(total_timer);
end


function error = relative_error(value,reference)
    error = norm(value-reference,'fro')/max(norm(reference,'fro'),realmin);
end
