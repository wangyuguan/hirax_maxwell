function [corrections,info] = fourfold_quad_corr_mats( ...
    S1,S2,S3,S4,eps_quad,zk,batch_size)
% FOURFOLD_QUAD_CORR_MATS Construct one block row of the C4 correction.
%
%  Syntax
%   [C,info] = fourfold_quad_corr_mats(S1,S2,S3,S4,eps_quad,zk,batch_size)
%
%  Input arguments:
%    * S1,S2,S3,S4: surfer objects in merge order, with matching node order
%    * eps_quad: requested quadrature accuracy
%    * zk: Helmholtz wavenumber
%    * batch_size: number of target rows per quadrature batch (2000)
%
%  Output arguments:
%    * C: sparse slp, gx, gy, gz blocks targeting S1, with C4 indexing data
%    * info: symmetry errors, block timings, and sparse storage counts
%
%  The clockwise rotation cycle is S1 -> S2 -> S4 -> S3 -> S1. Only the
%  geometry is assumed symmetric; the four densities may be unrelated.
%  Use fourfold_apply_quad_corr to apply the stored row to all four leaves.

    if nargin < 7 || isempty(batch_size)
      batch_size = 2000;
    end
    surfaces = {S1,S2,S4,S3};
    cycle_to_merge = [1,2,4,3];
    Q = [0,1,0;-1,0,0;0,0,1];

    npts = S1.npts;
    for isurface = 2:4
      assert(surfaces{isurface}.npts == npts, ...
          'All four surfaces must have the same number of nodes.')
      assert(isequal(surfaces{isurface}.norders,S1.norders) && ...
          isequal(surfaces{isurface}.iptype,S1.iptype), ...
          'Rotational reuse requires matching patch orders and types.')
    end

    scale = max(1,max(vecnorm(S1.r,2,1)));
    symmetry_tolerance = 1e3*eps(scale);
    position_errors = zeros(1,4);
    normal_errors = zeros(1,4);
    weight_errors = zeros(1,4);
    for isurface = 1:4
      inext = mod(isurface,4)+1;
      position_errors(isurface) = max(vecnorm( ...
          Q*surfaces{isurface}.r-surfaces{inext}.r,2,1));
      normal_errors(isurface) = max(vecnorm( ...
          Q*surfaces{isurface}.n-surfaces{inext}.n,2,1));
      weight_errors(isurface) = max(abs( ...
          surfaces{isurface}.wts-surfaces{inext}.wts));
    end
    assert(max(position_errors) < symmetry_tolerance && ...
        max(normal_errors) < symmetry_tolerance && ...
        max(weight_errors) < symmetry_tolerance, ...
        'The four surfaces do not have the required rotational symmetry.')

    % Construct one target block row. The self block uses singular
    % self-to-self quadrature; the other blocks are close-surface corrections.
    Bslp = cell(1,4);
    Bx = cell(1,4);
    By = cell(1,4);
    Bz = cell(1,4);
    block_seconds = zeros(1,4);
    for source_cycle = 1:4
      block_timer = tic;
      source_surface = surfaces{source_cycle};
      target = S1;
      if source_cycle == 1
        target = [];
      end
      B = nrccie_quad_corr_block(source_surface,eps_quad,zk,target,batch_size);
      [Bslp{source_cycle},Bx{source_cycle},By{source_cycle},Bz{source_cycle}] = B{:};
      block_seconds(source_cycle) = toc(block_timer);
    end

    % Surfaces on opposite corners have no near interaction at all, so the
    % matching block row entry is structurally empty and its multiply is
    % skipped when the correction is applied.
    nonempty = false(1,4);
    for source_cycle = 1:4
      nonempty(source_cycle) = nnz(Bslp{source_cycle})+ ...
          nnz(Bx{source_cycle})+nnz(By{source_cycle})+ ...
          nnz(Bz{source_cycle}) > 0;
    end

    corrections.slp = Bslp;
    corrections.gx = Bx;
    corrections.gy = By;
    corrections.gz = Bz;
    corrections.nonempty = nonempty;
    corrections.cycle_to_merge = cycle_to_merge;
    corrections.rotation = Q;
    corrections.component_npts = npts;
    corrections.npts = 4*npts;

    info.cycle_to_merge = cycle_to_merge;
    info.rotation = Q;
    info.position_errors = position_errors;
    info.normal_errors = normal_errors;
    info.weight_errors = weight_errors;
    info.nonempty = nonempty;
    info.block_seconds = block_seconds; % self, relative offsets 1, 2, 3
    info.base_block_nonzeros = zeros(4,4);
    for source_cycle = 1:4
      info.base_block_nonzeros(source_cycle,:) = [ ...
          nnz(Bslp{source_cycle}),nnz(Bx{source_cycle}), ...
          nnz(By{source_cycle}),nnz(Bz{source_cycle})];
    end
    info.stored_nonzeros = sum(info.base_block_nonzeros,'all');
    info.assembled_nonzeros = 4*info.stored_nonzeros;
end
