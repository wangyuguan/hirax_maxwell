function [array,errors] = toy_array_coupling( ...
    outer,model,Sref,Aref,Cref,incidentfun,targets,opts)
%TOY_ARRAY_COUPLING Reuse one station scattering matrix in a four-body array.
%   [array,errors] = toy_array_coupling(outer,model,Sref,Aref,Cref, ...
%       incidentfun,targets,opts)
%
%   outer contains Scat, proxy, input_projection, and output_interpolation.
%   model.station_shifts contains four translations of the same station.
%   incidentfun(r) returns interleaved [Ex,Ey,Ez,Hx,Hy,Hz] samples, with one
%   incident field per column. At most opts.array_nrhs columns are tested.
%
%   The compressed array solves b_i-S*sum_{j~=i} T_ij*b_j=S*a_i, where
%   T_ij=U'*L_i*Q_ji*U. Every station uses exactly the same square S. Each
%   Q_ji is evaluated in target batches, never as a full surface matrix.
%
%   The independent reference is the true corrected full-array NRCCIE,
%   applied by FMM and solved by GMRES with an exact single-station block
%   preconditioner. It retains the charge-continuity cross terms; a field-
%   only -R*K matrix is not used as the reference. Cref is the single-
%   station cell {Cslp,Cx,Cy,Cz}; its four translated copies are reused.
%   No ground response is included. Accuracy is relative to this physical
%   discretization, not a claim of mesh convergence.

    batch_size = option(opts,'array_batch_size',64);
    requested_rhs = option(opts,'array_nrhs',4);
    gmres_tol = option(opts,'array_gmres_tol',max(1e-9,20*opts.eps_fmm));
    gmres_restart = option(opts,'array_gmres_restart',60);
    gmres_cycles = option(opts,'array_gmres_cycles',8);
    accuracy_limit = option(opts,'array_accuracy_limit',5e-3);
    shifts = model.station_shifts;
    assert(isequal(size(shifts),[3,4]));
    assert(iscell(Cref) && numel(Cref)==4);
    assert(isequal(size(Aref),[3*Sref.npts,3*Sref.npts]));
    P = outer.proxy;
    U = outer.output_interpolation;
    V = outer.input_projection;
    Scat = outer.Scat;
    m = size(Scat,1);
    assert(isequal(size(Scat),[m,m]));
    assert(isequal(size(U),[4*P.npts,m]));
    assert(isequal(size(V),[m,4*P.npts]));
    timer = tic;
    timings = struct();

    % Translation invariance reduces twelve directed interactions to eight
    % distinct offsets for this square array. Construct each map only once.
    pairs = zeros(12,2);
    offsets = zeros(12,3);
    count = 0;
    for i = 1:4
      for j = 1:4
        if i==j, continue; end
        count = count+1;
        pairs(count,:) = [i,j];
        offsets(count,:) = (shifts(:,i)-shifts(:,j)).';
      end
    end
    [unique_offsets,~,offset_id] = unique(offsets,'rows');
    proxy = cell(1,4);
    for i = 1:4
      proxy{i} = translated_proxy(P,shifts(:,i));
    end
    EHi = incidentfun(proxy{1}.r);
    nrhs = min(requested_rhs,size(EHi,2));
    assert(nrhs>=1);
    rhs = complex(zeros(4*m,nrhs));
    external = rhs;
    for i = 1:4
      if i~=1, EHi = incidentfun(proxy{i}.r); end
      rows = (i-1)*m+(1:m);
      external(rows,:) = V*(maxwell_trace_matrix(proxy{i})*EHi(:,1:nrhs));
      rhs(rows,:) = Scat*external(rows,:);
    end
    clear EHi
    stage = tic;
    H = eye(4*m,'like',Scat);
    for d = 1:size(unique_offsets,1)
      T = translated_proxy(P,unique_offsets(d,:).');
      translation = compressed_translation(P,T,U,V,opts.zk,batch_size);
      block = -Scat*translation;
      for pair = find(offset_id==d).'
        i = pairs(pair,1); j = pairs(pair,2);
        H((i-1)*m+(1:m),(j-1)*m+(1:m)) = block;
      end
      fprintf('Array proxy translation %d/%d complete\n',d,size(unique_offsets,1));
    end
    clear translation block
    timings.translation = toc(stage);
    stage = tic;
    % The identity diagonal is already a natural block preconditioner.
    % Avoid a cubic factorization of the full compressed array matrix.
    b = complex(zeros(size(rhs)));
    compressed_flags = zeros(1,nrhs);
    compressed_iterations = zeros(nrhs,2);
    for j = 1:nrhs
      [b(:,j),compressed_flags(j),~,compressed_iterations(j,:)] = ...
          gmres(H,rhs(:,j),min(80,4*m),1e-12,8,[],[],rhs(:,j));
    end
    errors.compressed_solve = column_errors(H*b,rhs);
    timings.compressed_solve = toc(stage);
    clear H

    % Keep the reference independent of the compressed interface. Every
    % off-diagonal physical block is handled by the actual NRCCIE apply.
    stage = tic;
    stations = cell(1,4);
    for i = 1:4
      stations{i} = affine_transf(Sref,eye(3),shifts(:,i));
    end
    Sarray = merge([stations{:}]);
    parts = cell(1,4);
    for kernel = 1:4
      parts{kernel} = cell(4,4);
      for i = 1:4
        parts{kernel}{i,i} = Cref{kernel};
      end
    end
    cross_nonzeros = zeros(size(unique_offsets,1),4);
    for d = 1:size(unique_offsets,1)
      T = affine_transf(Sref,eye(3),unique_offsets(d,:).');
      C = nrccie_quad_corr_block(Sref,opts.eps_quad,opts.zk,T,batch_size);
      for kernel = 1:4
        cross_nonzeros(d,kernel) = nnz(C{kernel});
        for pair = find(offset_id==d).'
          i = pairs(pair,1); j = pairs(pair,2);
          parts{kernel}{i,j} = C{kernel};
        end
      end
      fprintf('Array reference correction %d/%d complete\n',d,size(unique_offsets,1));
    end
    Carray = cell(1,4);
    for kernel = 1:4
      Carray{kernel} = cell2mat(parts{kernel});
    end
    clear parts C T stations
    timings.reference_quadrature = toc(stage);
    Afun = @(q) nrccie_apply(Sarray,q,Carray,opts);
    Fstation = decomposition(Aref,'lu');
    station_dofs = size(Aref,1);
    Mfun = @(x) station_preconditioner(Fstation,x,station_dofs);
    EH = incidentfun(Sarray.r);
    physical_rhs = nrccie_rhs_matrix(Sarray,opts.alpha)*EH(:,1:nrhs);
    clear EH
    qref = complex(zeros(3*Sarray.npts,nrhs));
    flags = zeros(1,nrhs);
    iterations = zeros(nrhs,2);
    reported_residual = zeros(1,nrhs);
    residual_histories = cell(1,nrhs);
    errors.physical_reference_residual = zeros(1,nrhs);
    stage = tic;
    for j = 1:nrhs
      [qref(:,j),flags(j),reported_residual(j),iterations(j,:), ...
          residual_histories{j}] = gmres(Afun,physical_rhs(:,j), ...
          gmres_restart,gmres_tol,gmres_cycles,Mfun);
      errors.physical_reference_residual(j) = ...
          column_errors(Afun(qref(:,j)),physical_rhs(:,j));
      fprintf('Full-array NRCCIE RHS %d: flag %d, true residual %.3e\n', ...
          j,flags(j),errors.physical_reference_residual(j));
    end
    timings.reference_solve = toc(stage);

    stage = tic;
    ntargets = size(targets.r,2);
    field = complex(zeros(6*ntargets,nrhs));
    field_no_coupling = field;
    field_ref = field;
    for first = 1:batch_size:ntargets
      ids = first:min(first+batch_size-1,ntargets);
      rows = reshape(6*(ids-1)+(1:6).',1,[]);
      T = struct('r',targets.r(:,ids));
      for i = 1:4
        dofs = (i-1)*m+(1:m);
        Q = stratton_chu_matrix(proxy{i},T,opts.zk);
        field(rows,:) = field(rows,:)+Q*(U*b(dofs,:));
        field_no_coupling(rows,:) = field_no_coupling(rows,:)+Q*(U*rhs(dofs,:));
      end
      field_ref(rows,:) = nrccie_matrix(Sarray,T,opts.zk)*qref;
    end
    errors.target_vs_physical_reference = column_errors(field,field_ref);
    errors.no_array_coupling_vs_physical_reference = ...
        column_errors(field_no_coupling,field_ref);
    for kind = 1:2
      rows = reshape(6*(0:ntargets-1)+(3*(kind-1)+(1:3)).',1,[]);
      if kind==1
        errors.electric_target = column_errors(field(rows,:),field_ref(rows,:));
      else
        errors.magnetic_target = column_errors(field(rows,:),field_ref(rows,:));
      end
    end
    timings.field_validation = toc(stage);
    timings.total = toc(timer);
    stats = struct('station_scattering_constructions',1,'translated_reuses',3, ...
        'station_proxy_dofs',m,'compressed_array_dofs',4*m, ...
        'physical_array_dofs',3*Sarray.npts,'physical_array_nodes',Sarray.npts, ...
        'directed_interactions',12,'unique_translation_constructions',size(unique_offsets,1), ...
        'reference_cross_correction_constructions',size(unique_offsets,1), ...
        'reference_cross_correction_nonzeros',cross_nonzeros,'nrhs',nrhs, ...
        'ground_response_included',false);
    reference = struct('gmres_flags',flags,'gmres_iterations',iterations, ...
        'gmres_reported_residual',reported_residual,'gmres_tolerance',gmres_tol, ...
        'gmres_residual_histories',{residual_histories}, ...
        'description','Full corrected physical NRCCIE with all charge-continuity cross terms');
    array = struct('translated_proxies',{proxy},'station_shifts',shifts, ...
        'outgoing_coefficients',b,'external_incoming_coefficients',external, ...
        'compressed_gmres_flags',compressed_flags, ...
        'compressed_gmres_iterations',compressed_iterations, ...
        'field',field,'field_ref',field_ref, ...
        'field_no_coupling',field_no_coupling,'targets',targets, ...
        'stats',stats,'timings',timings,'reference',reference);
    array.accuracy_limit = accuracy_limit;
    array.passed = all(flags==0) && all(compressed_flags==0) && all(isfinite(field(:))) ...
        && max(errors.physical_reference_residual)<=1e-8 ...
        && max(errors.compressed_solve)<=1e-10 ...
        && max(errors.target_vs_physical_reference)<=accuracy_limit;
    fprintf('Array validation: %d compressed / %d physical DOFs; field error %.3e\n', ...
        4*m,3*Sarray.npts,max(errors.target_vs_physical_reference));
end

function translation = compressed_translation(P,T,U,V,zk,batch_size)
    translation = complex(zeros(size(V,1),size(U,2)));
    for first = 1:batch_size:T.npts
      ids = first:min(first+batch_size-1,T.npts);
      rows = reshape(4*(ids-1)+(1:4).',1,[]);
      subset = proxy_subset(T,ids);
      trace = maxwell_trace_matrix(subset)*stratton_chu_matrix(P,subset,zk);
      translation = translation+V(:,rows)*(trace*U);
    end
end

function value = station_preconditioner(F,x,n)
    value = zeros(size(x),'like',x);
    for i = 1:4
      rows = (i-1)*n+(1:n);
      value(rows,:) = F\x(rows,:);
    end
end

function T = translated_proxy(P,shift)
    T = P;
    T.r = P.r+shift;
end

function P = proxy_subset(P,ids)
    for name = {'r','n','dru','drv'}
      P.(name{1}) = P.(name{1})(:,ids);
    end
    P.wts = P.wts(ids);
    P.npts = numel(ids);
end

function value = option(opts,name,fallback)
    if isfield(opts,name) && ~isempty(opts.(name))
      value = opts.(name);
    else
      value = fallback;
    end
end

function error = column_errors(value,reference)
    error = vecnorm(value-reference,2,1)./max(vecnorm(reference,2,1),realmin);
end
