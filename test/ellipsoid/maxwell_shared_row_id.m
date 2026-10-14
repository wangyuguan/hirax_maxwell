function [nodes,T,info] = maxwell_shared_row_id(F,tol)
% MAXWELL_SHARED_ROW_ID Interpolate Maxwell traces from shared surface nodes.
%
%  Input arguments:
%    * F: 4*npts by nfields matrix of interleaved [Eu,Ev,Hu,Hv] traces
%        Include both incoming and outgoing fields, normalized by the caller.
%    * tol: relative Frobenius tolerance for the training matrix
%
%  Output arguments:
%    * nodes: selected surface node indices; all four components are retained
%    * T: interpolation matrix, F approximately equals T*F(info.dofs,:)
%        T(info.dofs,:) is the identity. Components may mix in interpolation.
%    * info: selected DOFs, numerical rank, residual, and coefficient size

    if nargin < 2, tol = 1e-6; end
    assert(mod(size(F,1),4)==0,'F must have four rows per surface node.');
    assert(isfinite(tol) && tol>0 && tol<1,'tol must lie between zero and one.');
    assert(all(isfinite(F(:))),'F must contain finite values.');
    npts = size(F,1)/4;
    G = F';
    residual = G;
    scale = norm(G,'fro');
    Q = zeros(size(G,1),min(size(G)),'like',G);
    available = true(1,npts);
    nodes = zeros(1,npts);
    history = zeros(1,npts+1);
    history(1) = double(scale>0);
    nselected = 0;
    numerical_rank = 0;

    % Block-pivoted QR of F': every pivot selects a complete surface node.
    % A small margin accommodates rounding in the final least-squares fit.
    while scale>0 && norm(residual,'fro')>0.5*tol*scale && any(available)
      scores = sum(reshape(sum(abs(residual).^2,1),4,npts),1);
      scores(~available) = -inf;
      [~,node] = max(scores);
      available(node) = false;
      block = residual(:,4*node-3:4*node);
      if numerical_rank>0
        basis = Q(:,1:numerical_rank);
        block = block-basis*(basis'*block);
        block = block-basis*(basis'*block);
      end
      [new_basis,R,~] = qr(block,0);
      diagonal = abs(diag(R));
      if isempty(diagonal) || max(diagonal)==0, continue; end
      threshold = max(size(block))*eps(class(F))*max(diagonal);
      rank_block = sum(diagonal>threshold);
      rank_block = min(rank_block,size(Q,2)-numerical_rank);
      if rank_block==0, continue; end
      new_basis = new_basis(:,1:rank_block);
      residual = residual-new_basis*(new_basis'*residual);
      Q(:,numerical_rank+(1:rank_block)) = new_basis;
      numerical_rank = numerical_rank+rank_block;
      nselected = nselected+1;
      nodes(nselected) = node;
      history(nselected+1) = norm(residual,'fro')/scale;
      if mod(nselected,50)==0
        fprintf('Joint ID: %d points, relative residual %.3e\n', ...
            nselected,history(nselected+1));
      end
      if numerical_rank==size(Q,2), break; end
    end

    nodes = nodes(1:nselected);
    dofs = reshape(4*(nodes-1)+(1:4).',1,[]);
    if nselected==0
      T = zeros(size(F,1),0,'like',F);
      relative_residual = 0;
      maxabsT = 0;
    else
      % A selected block can be rank deficient. The minimum-norm fit handles
      % dependent components without discarding any component at that node.
      T = lsqminnorm(G(:,dofs),G)';
      T(dofs,:) = eye(4*nselected,'like',F);
      relative_residual = norm(F-T*F(dofs,:),'fro')/scale;
      maxabsT = max(abs(T(:)));
    end
    info = struct('dofs',dofs,'npoints',nselected,'rank',numerical_rank, ...
        'relative_residual',relative_residual,'maxabsT',maxabsT, ...
        'residual_history',history(1:nselected+1));
end
