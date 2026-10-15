function results = check_toy_physics(norders,nproxies)
%CHECK_TOY_PHYSICS Separate physical refinement from proxy compression.
%   results = check_toy_physics([8,10],[28,40])
%
%   Solve the full corrected five-body NRCCIE for two independent plane-
%   wave polarizations, with a self-block preconditioned FMM/GMRES solve.
%   Compare its direct scattered field with the Stratton-Chu representation
%   of the same field sampled on the enclosing ellipsoid. No ID or proxy
%   scattering matrix is constructed. Comparing physical orders measures
%   field change, while comparing proxy orders isolates surface quadrature.
%   Geometry, frequency, and target points match the toy array rehearsal.
%   The experiment is in free space; the ground is only a height reference.
%
%   Results are saved after each physical order in
%   test/results/toy_physics_refinement.mat. Refinement differences are
%   evidence of convergence, not a bound against an exact solution.

    if nargin<1 || isempty(norders), norders = [8,10]; end
    if nargin<2 || isempty(nproxies), nproxies = [28,40]; end
    validateattributes(norders,{'numeric'},{'vector','integer','>=',4});
    validateattributes(nproxies,{'numeric'},{'vector','integer','>=',8});
    test_dir = fileparts(mfilename('fullpath'));
    run(fullfile(test_dir,'..','..','startup.m'));
    output_dir = fullfile(test_dir,'..','results');
    if ~isfolder(output_dir), mkdir(output_dir); end
    output_file = fullfile(output_dir,'toy_physics_refinement.mat');
    results = struct([]);
    for order = norders(:).'
      timer = tic;
      [model,physical] = toy_antenna_array_geometry(order);
      S = physical.station;
      opts = struct('zk',model.zk,'alpha',1,'jump',0.5,'eps_fmm',1e-11);
      C = nrccie_quad_corr_block(S,1e-10,model.zk,[],128);
      ncomponent = physical.body.npts;
      nd = 3*ncomponent;
      bodies = {physical.body,physical.leaves{1}};
      factors = cell(1,2);
      for j = 1:2
        ids = (j-1)*ncomponent+(1:ncomponent);
        localC = cellfun(@(c)c(ids,ids),C,'UniformOutput',false);
        A = nrccie_system_block(bodies{j},[],localC,opts,1:nd,1:nd);
        factors{j} = decomposition(A,'lu');
      end
      clear A
      Afun = @(q) nrccie_apply(S,q,C,opts);
      Mfun = @(q) prec(factors,q,nd);
      rhs = nrccie_rhs_matrix(S,1)*waves(S.r,model.zk);
      q = complex(zeros(size(rhs)));
      flags = zeros(1,2); true_residual = zeros(1,2);
      for wave = 1:2
        [q(:,wave),flags(wave)] = gmres(Afun,rhs(:,wave),60,1e-10,8,Mfun);
        true_residual(wave) = norm(Afun(q(:,wave))-rhs(:,wave))/norm(rhs(:,wave));
        fprintf('DIAG p%d wave%d flag%d residual%.4e time%.1f\n',order,wave,flags(wave),true_residual(wave),toc(timer));
      end
      direction = directions(36,0.37);
      targets.r = model.outer_center+[1.2*direction,1.8*direction,2.5*direction];
      direct_field = nrccie_matrix(S,targets,model.zk)*q;
      proxy_orders = nproxies(:).';
      green_errors = zeros(2,numel(proxy_orders));
      for j = 1:numel(proxy_orders)
        P = maxwell_ellipsoid(model.outer_axes,proxy_orders(j),model.outer_center);
        trace = complex(zeros(4*P.npts,2));
        for first = 1:64:P.npts
          ids = first:min(first+63,P.npts);
          T = subset(P,ids);
          rows = reshape(4*(ids-1)+(1:4).',[],1);
          trace(rows,:) = maxwell_trace_matrix(T)*(nrccie_matrix(S,T,model.zk)*q);
        end
        field_sc = stratton_chu_matrix(P,targets,model.zk)*trace;
        green_errors(:,j) = (vecnorm(field_sc-direct_field)./vecnorm(direct_field)).';
        fprintf('DIAG p%d q%d GREEN %.6e %.6e\n',order,proxy_orders(j),green_errors(:,j));
      end
      index = numel(results)+1;
      results(index).order = order;
      results(index).proxy_orders = proxy_orders;
      results(index).green_errors = green_errors;
      results(index).true_residual = true_residual;
      results(index).flags = flags;
      results(index).field = direct_field;
      results(index).seconds = toc(timer);
      if index>1
        results(index).field_change = vecnorm(direct_field-results(index-1).field)./vecnorm(direct_field);
        fprintf('DIAG p%d vs p%d FIELD CHANGE %.6e %.6e\n',order,results(index-1).order,results(index).field_change);
      end
      save(output_file,'results');
      fprintf('DIAG p%d SAVED seconds%.1f\n',order,toc(timer));
      clear physical S Afun Mfun C factors q rhs trace field_sc P
    end
    assert(all(arrayfun(@(r)all(r.flags==0) && max(r.true_residual)<1e-8,results)), ...
        'Physical reference solves did not reach the requested accuracy.');
    fprintf('TOY_PHYSICS_DIAGNOSTIC_DONE: %s\n',output_file);
end

function x = prec(factors,b,n)
    x = zeros(size(b),'like',b);
    x(1:n,:) = factors{1}\b(1:n,:);
    for j = 2:5
      ids = (j-1)*n+(1:n);
      x(ids,:) = factors{2}\b(ids,:);
    end
end

function EH = waves(r,k)
    d = directions(4,0.73);
    d = d(:,1);
    [~,idx] = min(abs(d));
    e = zeros(3,1); e(idx)=1;
    u = cross(d,e); u=u/norm(u); v=cross(d,u);
    phase = exp(1i*k*(d.'*r));
    EH = [reshape([u*phase;v*phase],[],1),reshape([v*phase;-u*phase],[],1)];
end

function d = directions(n,angle)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5))+angle;
    d = [sqrt(1-z.^2).*cos(phi);sqrt(1-z.^2).*sin(phi);z];
end

function P = subset(P,ids)
    for field = {'r','n','dru','drv'}
      P.(field{1}) = P.(field{1})(:,ids);
    end
    P.wts=P.wts(ids); P.npts=numel(ids);
end
