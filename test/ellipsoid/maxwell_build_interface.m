function [interface,BI,CI] = maxwell_build_interface(P,mapfun,A,tol,leaf_points,label,solvefun)
% MAXWELL_BUILD_INTERFACE Construct a compressed square scattering matrix.
%
%  Input arguments:
%    * P: closed proxy surface with r, n, dru, drv, wts, and npts
%    * mapfun: [B,C]=mapfun(ids) returns incoming columns B(:,dofs) and
%        outgoing rows C(dofs,:) for the requested proxy point indices.
%        There are four consecutive Eu,Ev,Hu,Hv components per point.
%    * A: square physical or hybrid system matrix, used to check residuals
%    * tol: relative joint-ID tolerance (default 1e-6)
%    * leaf_points: maximum points in a spatial ID leaf (default 64)
%    * label: progress-output label (default 'Proxy interface')
%    * solvefun: optional X=solvefun(rhs) applying A^{-1}. This may exploit
%        symmetry. If omitted, one LU decomposition of A is reused.
%
%  Output arguments:
%    * interface: proxy, shared skeleton nodes, input_projection=U',
%        output_interpolation=U, Scat, ID statistics, and measured errors
%    * BI,CI: optional selected maps, useful for validation or recovery
%
%  A deterministic multilevel ID of [B'/normB,C/normC] gives a common U:
%    B approximately equals BI*U', C approximately equals U*CI.
%  Only Scat=CI*(A\BI) is formed, of size 4*nskeleton by 4*nskeleton.
%  A new incident trace din is compressed as U'*din, not din at the nodes;
%  its outgoing trace is U*Scat*(U'*din). The full proxy scattering matrix
%  C*(A\B) is never assembled. Candidate B/C blocks and final error checks
%  are streamed; U is explicit in this small-problem implementation.
%
%  The returned object does not retain A, its factorization, BI, or CI.
%  errors.selected_residual measures A*X-BI for the selected solves.
%  errors.b and errors.c measure the complete maps in Frobenius norm.

    if nargin<4 || isempty(tol), tol = 1e-6; end
    if nargin<5 || isempty(leaf_points), leaf_points = 64; end
    if nargin<6 || isempty(label), label = 'Proxy interface'; end
    if nargin<7, solvefun = []; end
    validateattributes(tol,{'numeric'},{'scalar','real','finite','positive','<',1});
    validateattributes(leaf_points,{'numeric'},{'scalar','integer','positive'});
    assert(isa(mapfun,'function_handle'),'mapfun must be a function handle.');
    assert(isnumeric(A) && ismatrix(A) && size(A,1)==size(A,2) ...
        && ~isempty(A) && all(isfinite(A(:))),'A must be a finite square matrix.');
    assert(isempty(solvefun) || isa(solvefun,'function_handle'), ...
        'solvefun must be empty or a function handle.');
    timer_total = tic;
    timer_stage = tic;
    bnorm2 = 0;
    cnorm2 = 0;
    for first = 1:leaf_points:P.npts
      [B,C] = checked_maps(first:min(first+leaf_points-1,P.npts));
      bnorm2 = bnorm2+norm(B,'fro')^2;
      cnorm2 = cnorm2+norm(C,'fro')^2;
    end
    bnorm = sqrt(bnorm2);
    cnorm = sqrt(cnorm2);
    assert(bnorm>0 && cnorm>0,'Incoming and outgoing maps must be nonzero.');
    timings.normalization = toc(timer_stage);
    fprintf('%s: streamed norms; starting multilevel joint ID\n',label);

    timer_stage = tic;
    [nodes,U,info] = maxwell_multilevel_row_id(P,@joint_rows,tol,leaf_points);
    [BI,CI] = checked_maps(nodes);
    timings.id = toc(timer_stage);
    timer_stage = tic;
    berror2 = 0;
    cerror2 = 0;
    for first = 1:leaf_points:P.npts
      ids = first:min(first+leaf_points-1,P.npts);
      rows = point_dofs(ids);
      [B,C] = checked_maps(ids);
      berror2 = berror2+norm(B-BI*U(rows,:)','fro')^2;
      cerror2 = cerror2+norm(C-U(rows,:)*CI,'fro')^2;
    end
    errors = struct('b',sqrt(berror2/bnorm2),'c',sqrt(cerror2/cnorm2));
    dofs = point_dofs(nodes);
    assert(norm(U(dofs,:)-eye(numel(dofs)),'fro')<1e-11, ...
        'The selected trace block must be the identity.');
    timings.operator_check = toc(timer_stage);

    timer_stage = tic;
    if isempty(solvefun)
      F = decomposition(A,'lu');
      solvefun = @(rhs) F\rhs;
    end
    Scat = complex(zeros(size(CI,1),size(BI,2)));
    residual2 = 0;
    for first = 1:128:size(BI,2)
      ids = first:min(first+127,size(BI,2));
      rhs = BI(:,ids);
      X = solvefun(rhs);
      assert(isequal(size(X),size(rhs)) && all(isfinite(X(:))), ...
          'solvefun returned incompatible or nonfinite values.');
      Scat(:,ids) = CI*X;
      residual2 = residual2+norm(A*X-rhs,'fro')^2;
    end
    errors.selected_residual = sqrt(residual2)/norm(BI,'fro');
    timings.selected_solve = toc(timer_stage);
    timings.total = toc(timer_total);
    interface = struct('proxy',P,'nodes',nodes,'dofs',dofs, ...
        'skeleton',proxy_subset(P,nodes),'input_projection',U', ...
        'output_interpolation',U,'Scat',Scat,'id_info',info,'errors',errors, ...
        'normalization',struct('b_frobenius',bnorm,'c_frobenius',cnorm), ...
        'timings',timings, ...
        'trace_order','Eu,Ev,Hu,Hv per node; unit tangents; no weights');
    fprintf('%s: %d/%d points; B/C errors %.3e / %.3e; solve %.3e\n', ...
        label,numel(nodes),P.npts,errors.b,errors.c,errors.selected_residual);
    fprintf('%s: square scattering matrix %d by %d\n',label,size(Scat));

    function [B,C] = checked_maps(ids)
        [B,C] = mapfun(ids);
        assert(isequal(size(B),[size(A,1),4*numel(ids)]) ...
            && isequal(size(C),[4*numel(ids),size(A,2)]), ...
            'mapfun must return system-by-trace B and trace-by-system C.');
        assert(all(isfinite(B(:))) && all(isfinite(C(:))), ...
            'mapfun returned nonfinite values.');
    end

    function rows = joint_rows(ids)
        [B,C] = checked_maps(ids);
        rows = [B'/bnorm,C/cnorm];
    end
end

function dofs = point_dofs(nodes)
    dofs = reshape(4*(nodes(:).'-1)+(1:4).',1,[]);
end

function P = proxy_subset(P,ids)
    for name = {'r','n','dru','drv'}
      P.(name{1}) = P.(name{1})(:,ids);
    end
    P.wts = P.wts(ids);
    P.npts = numel(ids);
end
