function [A,info,cache] = toy_fourfold_system(leaves,opts,eps_quad,run_reference)
%TOY_FOURFOLD_SYSTEM Build and solve the toy feed using its C4 symmetry.
%   [A,info,cache] = toy_fourfold_system(leaves,opts,eps_quad,true)
%
%   leaves is the four-element cell array returned by
%   toy_antenna_array_geometry, ordered counterclockwise around the z axis.
%   Corresponding nodes and tangent frames must be rotations of leaf 1.
%   opts contains zk, alpha, and jump; optional batch_size defaults to 2000.
%
%   Only four quadrature blocks and four corrected NRCCIE blocks are built:
%   one self block and three relative interleaf blocks targeting leaf 1.
%   A is assembled in the original merge order for this small experiment.
%   cache.solve(rhs) uses four smaller Fourier-mode LU factors, not an LU
%   factorization of A. All right-hand sides are supported, without imposing
%   symmetry on the incident field or unknown densities.
%
%   cache.corrections is compatible with fourfold_apply_quad_corr in the
%   ORIGINAL leaf order. cache.system_blocks{r+1} is the block targeting
%   leaf 1 from leaf 1+r. If run_reference is true, independent quadrature
%   on the complete merged surface checks all sixteen reused blocks.

    if nargin<3 || isempty(eps_quad), eps_quad = 1e-10; end
    if nargin<4 || isempty(run_reference), run_reference = false; end
    assert(iscell(leaves) && numel(leaves)==4, ...
        'toy_fourfold_system:Leaves','Expected four leaves in CCW order.');
    assert(all(isfield(opts,{'zk','alpha','jump'})), ...
        'toy_fourfold_system:Options','opts requires zk, alpha, and jump.');
    batch_size = 2000;
    if isfield(opts,'batch_size'), batch_size = opts.batch_size; end
    npts = leaves{1}.npts;
    n = 3*npts;
    Q = [0,-1,0;1,0,0;0,0,1];
    frame_errors = zeros(4,2);
    for j = 1:4
      next = mod(j,4)+1;
      assert(leaves{j}.npts==npts, ...
          'toy_fourfold_system:Size','All leaves need the same node count.');
      frame_errors(j,1) = max(vecnorm(Q*leaves{j}.dru-leaves{next}.dru,2,1));
      frame_errors(j,2) = max(vecnorm(Q*leaves{j}.drv-leaves{next}.drv,2,1));
    end
    assert(max(frame_errors,[],'all')<1e-11, ...
        'toy_fourfold_system:Frames','Tangent frames do not rotate with C4.');

    % Existing leaf utilities use clockwise cycle [1,2,4,3]. Present our
    % CCW leaves as [1,4,2,3], then restore their original indices in C4.
    stage = tic;
    [C4,quad_info] = fourfold_quad_corr_mats( ...
        leaves{1},leaves{4},leaves{2},leaves{3},eps_quad,opts.zk,batch_size);
    C4.cycle_to_merge = [1,4,3,2];
    quad_info.cycle_to_merge = C4.cycle_to_merge;
    info.quadrature_seconds = toc(stage);
    info.quadrature = quad_info;
    info.tangent_frame_errors = frame_errors;
    info.quadrature_blocks_computed = 4;
    info.quadrature_blocks_without_reuse = 16;
    info.self_blocks_computed = 1;
    info.interleaf_blocks_computed = 3;
    info.matrix_blocks_computed = 4;
    info.matrix_blocks_without_reuse = 16;
    info.leaf_nodes = npts;
    info.leaf_dofs = n;
    info.feed_dofs = 4*n;

    stage = tic;
    blocks = cell(1,4);
    for source = 1:4
      relative = mod(1-source,4)+1; % CCW source -> clockwise correction
      C = {C4.slp{relative},C4.gx{relative}, ...
          C4.gy{relative},C4.gz{relative}};
      target = leaves{1};
      if source==1, target = []; end
      blocks{source} = nrccie_system_block( ...
          leaves{source},target,C,opts,1:n,1:n);
    end
    A = complex(zeros(4*n));
    for target = 1:4
      rows = (target-1)*n+(1:n);
      for source = 1:4
        cols = (source-1)*n+(1:n);
        A(rows,cols) = blocks{mod(source-target,4)+1};
      end
    end
    info.assembly_seconds = toc(stage);

    % A_ij=B_(j-i mod 4). With x_j=exp(2*pi*i*m*j/4)*xhat_m,
    % the mode matrix is sum_r exp(2*pi*i*m*r/4)*B_r.
    stage = tic;
    factors = cell(1,4);
    for mode = 0:3
      D = complex(zeros(n));
      for relative = 0:3
        D = D+exp(2*pi*1i*mode*relative/4)*blocks{relative+1};
      end
      factors{mode+1} = decomposition(D,'lu');
    end
    cache.solve = @(rhs) fourier_solve(factors,n,rhs);
    cache.corrections = C4;
    cache.system_blocks = blocks;
    info.factor_seconds = toc(stage);
    info.factor_sizes = n*ones(1,4);
    info.factor_entries_relative_to_full = 1/4;
    info.dense_factor_flops_relative_to_full = 1/16;

    % Deterministic, unrelated right-hand sides test all Fourier modes.
    indices = (1:4*n).';
    rhs = sin(indices*[0.19,0.37,0.71]) ...
        +1i*cos(indices*[0.23,0.43,0.83]);
    solution = cache.solve(rhs);
    info.fourier_solve_residual = norm(A*solution-rhs,'fro')/norm(rhs,'fro');
    assert(info.fourier_solve_residual<1e-10, ...
        'toy_fourfold_system:FourierSolve','Fourier factor solve failed.');

    info.reference_run = logical(run_reference);
    if run_reference
      stage = tic;
      S = merge([leaves{:}]);
      full_C = nrccie_quad_corr_block(S,eps_quad,opts.zk,[],batch_size);
      correction_error_sq = zeros(1,4);
      correction_norm_sq = zeros(1,4);
      matrix_error_sq = 0;
      matrix_norm_sq = 0;
      for target = 1:4
        node_rows = (target-1)*npts+(1:npts);
        rows = (target-1)*n+(1:n);
        % Stored derivatives target leaf 1; rotate to this target frame.
        Qt = Q^(target-1);
        for source = 1:4
          node_cols = (source-1)*npts+(1:npts);
          cols = (source-1)*n+(1:n);
          relative = mod(target-source,4)+1;
          reused = {C4.slp{relative}, ...
              Qt(1,1)*C4.gx{relative}+Qt(1,2)*C4.gy{relative}, ...
              Qt(2,1)*C4.gx{relative}+Qt(2,2)*C4.gy{relative}, ...
              C4.gz{relative}};
          for component = 1:4
            reference = full_C{component}(node_rows,node_cols);
            correction_error_sq(component) = correction_error_sq(component) ...
                +norm(reused{component}-reference,'fro')^2;
            correction_norm_sq(component) = correction_norm_sq(component) ...
                +norm(reference,'fro')^2;
          end
          reference_A = nrccie_system_block(S,[],full_C,opts,rows,cols);
          matrix_error_sq = matrix_error_sq+norm(A(rows,cols)-reference_A,'fro')^2;
          matrix_norm_sq = matrix_norm_sq+norm(reference_A,'fro')^2;
        end
      end
      info.reference_correction_errors = sqrt(correction_error_sq ...
          ./max(correction_norm_sq,realmin));
      info.reference_matrix_error = sqrt(matrix_error_sq/max(matrix_norm_sq,realmin));
      % Also exercise the existing compact correction apply independently
      % of the direct block-assembly route used to construct A above.
      node_indices = 1:4*npts;
      density = sin([0.17;0.29;0.53]*node_indices) ...
          +1i*cos([0.31;0.47;0.59]*node_indices);
      applied = cell(1,4);
      [applied{:}] = fourfold_apply_quad_corr(C4,density);
      info.reference_apply_errors = zeros(1,4);
      for component = 1:4
        reference = (full_C{component}*density.').';
        info.reference_apply_errors(component) = ...
            norm(applied{component}-reference,'fro')/max(norm(reference,'fro'),realmin);
      end
      info.reference_seconds = toc(stage);
      threshold = max(5e-8,100*eps_quad);
      assert(info.reference_matrix_error<threshold && ...
          max(info.reference_correction_errors)<threshold && ...
          max(info.reference_apply_errors)<threshold, ...
          'toy_fourfold_system:Reference','C4 reuse disagrees with full quadrature.');
    end
    fprintf('Toy feed C4: 4/16 correction blocks; four %d-by-%d mode factors\n',n,n);
    fprintf('C4 solve residual %.3e\n',info.fourier_solve_residual);
    if run_reference
      fprintf('Full-quadrature matrix comparison %.3e; correction max %.3e\n', ...
          info.reference_matrix_error,max(info.reference_correction_errors));
    end
end

function solution = fourier_solve(factors,n,rhs)
    assert(ismatrix(rhs) && size(rhs,1)==4*n, ...
        'toy_fourfold_system:RightHandSide','rhs must have 4*n rows.');
    nrhs = size(rhs,2);
    transformed = fft(reshape(rhs,n,4,nrhs),[],2);
    for mode = 1:4
      transformed(:,mode,:) = reshape( ...
          factors{mode}\reshape(transformed(:,mode,:),n,nrhs),n,1,nrhs);
    end
    solution = reshape(ifft(transformed,[],2),4*n,nrhs);
end
