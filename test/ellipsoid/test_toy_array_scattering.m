function result = test_toy_array_scattering(norder,nproxy,tol,run_array)
% TEST_TOY_ARRAY_SCATTERING Rehearse symmetry, nested compression, and reuse.
%   result = test_toy_array_scattering(8,28,1e-5,true)
%
%   Four rotated closed leaves -> one feed interface -> feed coupled to
%   the physical oblate body -> one whole-station ellipsoid interface ->
%   four translated stations. Ground height is preserved geometrically;
%   this experiment is in free space and includes no ground reflection.
%
%   Quadrature uses the existing compact C4 machinery. Four Fourier-mode
%   solves factor the feed exactly, including nonsymmetric incident fields.
%   Both interfaces use a deterministic multilevel common incoming/outgoing
%   point ID; no uncompressed proxy scattering matrix is formed.
%   Dense local matrices are intentional for this small validation problem.

    if nargin<1 || isempty(norder), norder = 8; end
    if nargin<2 || isempty(nproxy), nproxy = 28; end
    if nargin<3 || isempty(tol), tol = 1e-5; end
    if nargin<4 || isempty(run_array), run_array = true; end
    validateattributes(norder,{'numeric'},{'scalar','integer','>=',4});
    validateattributes(nproxy,{'numeric'},{'scalar','integer','>=',8});
    validateattributes(tol,{'numeric'},{'scalar','positive','<',1});
    total_timer = tic;
    [model,physical] = toy_antenna_array_geometry(norder);
    zk = model.zk;
    alpha = 1;
    opts = struct('zk',zk,'alpha',alpha,'jump',0.5,'eps_fmm',1e-10, ...
        'eps_quad',1e-10,'batch_size',128);
    Sf = physical.feed;
    Sd = physical.body;
    Sref = physical.station;
    nf = 3*Sf.npts;
    nd = 3*Sd.npts;
    Pf = maxwell_ellipsoid(model.inner_radius*ones(3,1),nproxy,model.inner_center);
    Po = maxwell_ellipsoid(model.outer_axes,nproxy,model.outer_center);
    Rf = nrccie_rhs_matrix(Sf,alpha);
    Rd = nrccie_rhs_matrix(Sd,alpha);
    Lf = maxwell_trace_matrix(Pf);
    Lo = maxwell_trace_matrix(Po);
    leaf_points = 64;
    incidentfun = @(r) incident_fields(r,zk,model.outer_center);

    fprintf('\nStage 1: four-leaf symmetry and feed scattering matrix\n');
    stage = tic;
    [Af,symmetry,feed_solver] = toy_fourfold_system( ...
        physical.leaves,opts,opts.eps_quad,false);
    timings.feed_system = toc(stage);
    feed_maps = @(ids) body_maps(ids,Pf,Sf,Rf,zk);
    stage = tic;
    [feed,BfI,CfI] = maxwell_build_interface(Pf,feed_maps,Af,tol, ...
        leaf_points,'Four-leaf sphere',feed_solver.solve);
    feed.zk = zk;
    feed.alpha = alpha;
    timings.feed_interface = toc(stage);
    mf = size(feed.Scat,1);

    fprintf('\nStage 2: couple the feed interface to the physical body\n');
    stage = tic;
    Cd = nrccie_quad_corr_block(Sd,opts.eps_quad,zk,[],opts.batch_size);
    Ad = nrccie_system_block(Sd,[],Cd,opts,1:nd,1:nd);
    Gfd = feed.input_projection*Lf*nrccie_matrix(Sd,Pf,zk);
    Gdf = Rd*stratton_chu_matrix(Pf,Sd,zk)*feed.output_interpolation;
    H = [Ad,-Gdf;-feed.Scat*Gfd,eye(mf)];
    FH = decomposition(H,'lu');
    timings.hybrid_system = toc(stage);
    feed_drive = feed.Scat*feed.input_projection*Lf;
    outer_maps = @(ids) hybrid_maps(ids,Po,Sd,Rd,Pf, ...
        feed.output_interpolation,feed_drive,zk);
    stage = tic;
    [outer,BoI,CoI] = maxwell_build_interface(Po,outer_maps,H,tol, ...
        leaf_points,'Whole-station ellipsoid',@(rhs) FH\rhs);
    outer.zk = zk;
    outer.alpha = alpha;
    timings.outer_interface = toc(stage);

    fprintf('\nStage 3: independent full five-body NRCCIE reference\n');
    stage = tic;
    % This full quadrature is for verification only, not interface construction.
    Cref = nrccie_quad_corr_block(Sref,opts.eps_quad,zk,[],opts.batch_size);
    Aref = nrccie_system_block(Sref,[],Cref,opts,1:nd+nf,1:nd+nf);
    Fref = decomposition(Aref,'lu');
    errors.symmetry_matrix = norm(Aref(nd+1:end,nd+1:end)-Af,'fro')/norm(Af,'fro');
    reference_feed = cellfun(@(C) C(Sd.npts+1:end,Sd.npts+1:end),Cref, ...
        'UniformOutput',false);
    density_test = reshape(sin(sqrt(2)*(1:4*Sf.npts))+ ...
        1i*cos(sqrt(3)*(1:4*Sf.npts)),4,Sf.npts);
    [p,gx,gy,gz] = fourfold_apply_quad_corr(feed_solver.corrections,density_test);
    applied = {p,gx,gy,gz};
    for j = 1:4
      expected = (reference_feed{j}*density_test.').';
      errors.symmetry_quadrature(j) = norm(applied{j}-expected,'fro')/max(norm(expected,'fro'),realmin);
    end
    [incd,labels] = incident_fields(Sd.r,zk,model.outer_center);
    incf = incidentfun(Sf.r);
    dinf = Lf*incidentfun(Pf.r);
    dino = Lo*incidentfun(Po.r);
    rhs = [Rd*incd;Rf*incf];
    qref = Fref\rhs;
    rhs_hybrid = [Rd*incd;feed.Scat*(feed.input_projection*dinf)];
    zh = FH\rhs_hybrid;
    qf_hybrid = feed_solver.solve(BfI*(feed.input_projection*dinf+Gfd*zh(1:nd,:)));
    qh = [zh(1:nd,:);qf_hybrid];
    errors.full_solve = relative_errors(Aref*qref,rhs);
    errors.hybrid_solve = relative_errors(H*zh,rhs_hybrid);
    errors.hybrid_full_residual = relative_errors(Aref*qh,rhs);
    errors.hybrid_density = relative_errors(qh,qref);
    errors.feed_recovery = relative_errors(CfI*qf_hybrid,zh(nd+1:end,:));
    errors.feed_incoming = relative_errors(BfI*(feed.input_projection*dinf),Rf*incf);
    errors.outer_incoming = relative_errors(BoI*(outer.input_projection*dino),rhs_hybrid);
    bouter = outer.Scat*(outer.input_projection*dino);
    errors.outer_trace_vs_hybrid = relative_errors( ...
        outer.output_interpolation*bouter,outer.output_interpolation*(CoI*zh));

    target_radii = [1.2,1.8,2.5];
    d = directions(36,0.37);
    targets.r = model.outer_center+[target_radii(1)*d,target_radii(2)*d,target_radii(3)*d];
    Kd = nrccie_matrix(Sd,targets,zk);
    Kf = nrccie_matrix(Sf,targets,zk);
    Qf = stratton_chu_matrix(Pf,targets,zk)*feed.output_interpolation;
    Qo = stratton_chu_matrix(Po,targets,zk);
    fields.reference = Kd*qref(1:nd,:)+Kf*qref(nd+1:end,:);
    fields.hybrid = Kd*zh(1:nd,:)+Qf*zh(nd+1:end,:);
    fields.outer = (Qo*outer.output_interpolation)*bouter;
    errors.hybrid_field = relative_errors(fields.hybrid,fields.reference);
    errors.outer_field = relative_errors(fields.outer,fields.reference);
    qf_alone = feed_solver.solve(Rf*incf);
    feed_b = feed.Scat*(feed.input_projection*dinf);
    errors.feed_field = relative_errors(Qf*feed_b,Kf*qf_alone);
    qd_once = Ad\(Rd*incd+Gdf*feed_b);
    fields.no_feedback = Kd*qd_once+Qf*feed_b;
    errors.no_feed_body_feedback = relative_errors(fields.no_feedback,fields.reference);

    % Independently assess Green representation, separate from both IDs.
    reference_trace = complex(zeros(4*Po.npts,numel(labels)));
    for first = 1:leaf_points:Po.npts
      ids = first:min(first+leaf_points-1,Po.npts);
      T = proxy_subset(Po,ids);
      rows = point_dofs(ids);
      reference_trace(rows,:) = maxwell_trace_matrix(T)*( ...
          nrccie_matrix(Sd,T,zk)*qref(1:nd,:)+nrccie_matrix(Sf,T,zk)*qref(nd+1:end,:));
    end
    errors.outer_green_reference = relative_errors(Qo*reference_trace,fields.reference);
    errors.outer_trace_vs_full = relative_errors(outer.output_interpolation*bouter,reference_trace);
    % NRCCIE has a charge-continuity term: compare it on solved densities.
    errors.feed_cross_compatibility = relative_errors( ...
        -Rd*nrccie_matrix(Sf,Sd,zk)*qref(nd+1:end,:),Aref(1:nd,nd+1:end)*qref(nd+1:end,:));
    errors.body_cross_compatibility = relative_errors( ...
        -Rf*nrccie_matrix(Sd,Sf,zk)*qref(1:nd,:),Aref(nd+1:end,1:nd)*qref(1:nd,:));
    timings.station_validation = toc(stage);

    stats = struct('feed_physical_dofs',nf,'station_physical_dofs',nd+nf, ...
        'hybrid_dofs',nd+mf,'feed_proxy_points',Pf.npts, ...
        'feed_skeleton_points',numel(feed.nodes),'outer_proxy_points',Po.npts, ...
        'outer_skeleton_points',numel(outer.nodes),'feed_scattering_size',size(feed.Scat), ...
        'outer_scattering_size',size(outer.Scat),'station_interface_builds',1, ...
        'translated_reuses',0,'ground_response_included',false);
    array = [];
    output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if ~isfolder(output_dir), mkdir(output_dir); end
    tag = sprintf('toy_array_p%d_q%d_%s',norder,nproxy,strrep(sprintf('%.0e',tol),'-','m'));
    output_file = fullfile(output_dir,[tag,'.mat']);
    save(output_file,'model','feed','outer','symmetry','stats','errors','timings', ...
        'norder','nproxy','tol','labels','fields','targets','run_array','-v7.3');
    fprintf('Station checkpoint saved: %s\n',output_file);
    print_errors(errors);

    if run_array
      fprintf('\nStage 4: reuse one interface for four coupled stations\n');
      clear Af Ad Cd H FH Fref feed_solver feed_maps outer_maps BfI CfI BoI CoI
      clear Kd Kf Qf Qo Gfd Gdf feed_drive reference_feed applied density_test
      clear qref qh qf_hybrid qf_alone qd_once zh reference_trace
      stage = tic;
      array_targets.r = model.outer_center+[4*directions(48,0.21),6*directions(48,0.43)];
      [array,array_errors] = toy_array_coupling(outer,model,Sref,Aref,Cref, ...
          incidentfun,array_targets,opts);
      errors.array = array_errors;
      stats.translated_reuses = array.stats.translated_reuses;
      timings.array = toc(stage);
    end
    timings.total = toc(total_timer);
    accuracy_limit = 5e-3;
    passed = errors.symmetry_matrix<1e-9 && max(errors.symmetry_quadrature)<1e-9 ...
        && feed.errors.selected_residual<1e-10 && outer.errors.selected_residual<1e-10 ...
        && max([feed.errors.b,feed.errors.c,outer.errors.b,outer.errors.c])<10*tol ...
        && numel(feed.nodes)<Pf.npts && numel(outer.nodes)<Po.npts ...
        && max(errors.feed_field)<accuracy_limit ...
        && max(errors.feed_incoming)<accuracy_limit ...
        && max(errors.outer_incoming)<accuracy_limit ...
        && max(errors.hybrid_full_residual)<accuracy_limit ...
        && max(errors.outer_field)<accuracy_limit;
    if run_array
      passed = passed && array.passed;
    end
    result = struct('passed',passed,'model',model,'feed',feed,'outer',outer, ...
        'symmetry',symmetry,'stats',stats,'errors',errors,'timings',timings, ...
        'norder',norder,'nproxy',nproxy,'tol',tol,'run_array',run_array, ...
        'labels',{labels},'fields',fields,'targets',targets,'array',array, ...
        'accuracy_limit',accuracy_limit,'output_file',output_file, ...
        'reference','Same-discretization corrected NRCCIE; no mesh-convergence claim');
    save(output_file,'-struct','result','-v7.3');
    fprintf('Saved %s; total %.1f s\n',output_file,timings.total);
    assert(passed,'Toy rehearsal did not meet validation thresholds; inspect saved errors.');
    fprintf('TOY_ARRAY_SCATTERING_PASS\n');
end

function [B,C] = body_maps(ids,P,S,R,zk)
    T = proxy_subset(P,ids);
    B = -R*stratton_chu_matrix(T,S,zk);
    C = maxwell_trace_matrix(T)*nrccie_matrix(S,T,zk);
end

function [B,C] = hybrid_maps(ids,P,S,R,Pf,Uf,feed_drive,zk)
    T = proxy_subset(P,ids);
    B = [-R*stratton_chu_matrix(T,S,zk);-feed_drive*stratton_chu_matrix(T,Pf,zk)];
    L = maxwell_trace_matrix(T);
    C = [L*nrccie_matrix(S,T,zk),L*stratton_chu_matrix(Pf,T,zk)*Uf];
end

function P = proxy_subset(P,ids)
    for field = {'r','n','dru','drv'}
      P.(field{1}) = P.(field{1})(:,ids);
    end
    P.wts = P.wts(ids);
    P.npts = numel(ids);
end

function dofs = point_dofs(ids)
    dofs = reshape(4*(ids-1)+(1:4).',1,[]);
end

function [EH,labels] = incident_fields(r,zk,center)
    d = directions(4,0.73);
    EH = complex(zeros(6*size(r,2),10));
    labels = cell(1,10);
    for j = 1:4
      [~,axis_id] = min(abs(d(:,j)));
      e = zeros(3,1); e(axis_id) = 1;
      u = cross(d(:,j),e); u = u/norm(u);
      v = cross(d(:,j),u);
      phase = exp(1i*zk*(d(:,j).'*r));
      EH(:,2*j-1) = reshape([u*phase;v*phase],[],1);
      EH(:,2*j) = reshape([v*phase;-u*phase],[],1);
      labels{2*j-1} = sprintf('plane_%d_u',j);
      labels{2*j} = sprintf('plane_%d_v',j);
    end
    source.r = center+[0.2;-0.3;3];
    source.edips = [1;0.3i;0.2];
    source.hdips = [0.2;1;0.1i];
    [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'ed');
    EH(:,9) = reshape([E;H],[],1); labels{9} = 'electric_dipole';
    [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'hd');
    EH(:,10) = reshape([E;H],[],1); labels{10} = 'magnetic_dipole';
end

function d = directions(n,angle)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5))+angle;
    d = [sqrt(1-z.^2).*cos(phi);sqrt(1-z.^2).*sin(phi);z];
end

function e = relative_errors(value,reference)
    e = vecnorm(value-reference,2,1)./max(vecnorm(reference,2,1),realmin);
end

function print_errors(errors)
    names = fieldnames(errors);
    fprintf('\nMaximum relative errors over independent incident fields:\n');
    for j = 1:numel(names)
      value = errors.(names{j});
      fprintf('  %-36s %.6e\n',names{j},max(value(:)));
    end
end
