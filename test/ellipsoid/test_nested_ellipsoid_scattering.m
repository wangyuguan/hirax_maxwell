function result = test_nested_ellipsoid_scattering(norder,nproxy,tol)
%TEST_NESTED_ELLIPSOID_SCATTERING Couple a compressed particle to a physical BIE.
%   result = test_nested_ellipsoid_scattering(12,20,1e-6)
%
%   A small PEC ellipsoid is replaced by its sphere scattering matrix. A
%   larger, physically disjoint ellipsoid retains its NRCCIE unknowns. One
%   outer sphere encloses both bodies and the small sphere. Its scattering
%   matrix represents the coupled pair, including repeated reflections.
%
%   Both interfaces use deterministic multilevel joint B/C point IDs, with
%   kernel callbacks for candidate blocks. No full sphere scattering matrix
%   is formed. Physical BIE and hybrid matrices are dense for this experiment.
%   Incoming coefficients are U'*din; they are not din at skeleton nodes.
%
%   Validation includes a full two-body corrected NRCCIE reference. The
%   full reference retains the charge-continuity cross terms; -R*K is checked
%   separately and is not substituted for the true NRCCIE cross blocks.

    if nargin<1 || isempty(norder), norder = 12; end
    if nargin<2 || isempty(nproxy), nproxy = 20; end
    if nargin<3 || isempty(tol), tol = 1e-6; end
    validateattributes(norder,{'numeric'},{'scalar','integer','>=',4});
    validateattributes(nproxy,{'numeric'},{'scalar','integer','>=',8});
    validateattributes(tol,{'numeric'},{'scalar','positive','<',1});
    total_timer = tic;
    zk = 2*pi;
    alpha = 1;
    eps_quad = 1e-10;
    opts = struct('zk',zk,'alpha',alpha,'jump',0.5);
    geometry.small_axes = [0.12;0.10;0.08];
    geometry.small_center = [0;0;0];
    geometry.large_axes = [0.35;0.28;0.22];
    geometry.large_center = [0.90;0;0];
    geometry.small_radius = 0.30;
    geometry.outer_center = [0.45;0;0];
    geometry.outer_radius = 1.20;
    small_proxy_order = max(12,nproxy-4);
    leaf_points = 64;

    Ss = geometries.ellipsoid(geometry.small_axes,[1;1;1], ...
        geometry.small_center,norder,11);
    Sd = geometries.ellipsoid(geometry.large_axes,[1;1;1], ...
        geometry.large_center,norder,11);
    Ps = maxwell_sphere(geometry.small_radius,small_proxy_order);
    Ps.r = Ps.r+geometry.small_center;
    Po = maxwell_sphere(geometry.outer_radius,nproxy);
    Po.r = Po.r+geometry.outer_center;
    Rs = nrccie_rhs_matrix(Ss,alpha);
    Rd = nrccie_rhs_matrix(Sd,alpha);
    Ls = maxwell_trace_matrix(Ps);
    Lo = maxwell_trace_matrix(Po);
    ns = 3*Ss.npts;
    nd = 3*Sd.npts;
    assert(norm(geometry.large_center-geometry.small_center) ...
        > max(geometry.large_axes)+geometry.small_radius);
    assert(norm(geometry.small_center-geometry.outer_center) ...
        +geometry.small_radius<geometry.outer_radius);
    assert(norm(geometry.large_center-geometry.outer_center) ...
        +max(geometry.large_axes)<geometry.outer_radius);
    fprintf('Nested ellipsoids: %d small + %d large physical unknowns\n',ns,nd);
    fprintf('Sphere points: %d small / %d outer; physical order %d\n', ...
        Ps.npts,Po.npts,norder);

    stage = tic;
    As = physical_matrix(Ss,[],opts,eps_quad);
    Fs = decomposition(As,'lu');
    timings.small_physical = toc(stage);
    stage = tic;
    small_maps = @(ids) body_maps(ids,Ps,Ss,Rs,zk);
    small = build_interface(Ps,small_maps,tol,leaf_points,'Small sphere');
    [BsI,CsI] = body_maps(small.nodes,Ps,Ss,Rs,zk);
    [small.Scat,errors.small_selected_residual] = ...
        selected_response(As,Fs,BsI,CsI);
    small.zk = zk;
    small.alpha = alpha;
    timings.small_interface = toc(stage);
    ms = size(small.Scat,1);
    fprintf('Small scattering matrix: %d by %d\n',ms,ms);

    stage = tic;
    Ad = physical_matrix(Sd,[],opts,eps_quad);
    % The small component is now accessed only through Ps, U_s, and Scat_s.
    % These are regular interior and outgoing exterior sphere evaluations.
    large_to_small_trace = Ls*nrccie_matrix(Sd,Ps,zk);
    Gsd = small.input_projection*large_to_small_trace;
    Gds = Rd*stratton_chu_matrix(Ps,Sd,zk)*small.output_interpolation;
    H = [Ad,-Gds; -small.Scat*Gsd,eye(ms)];
    FH = decomposition(H,'lu');
    timings.hybrid_system = toc(stage);
    fprintf('Hybrid system: %d physical large + %d small interface DOFs\n',nd,ms);

    stage = tic;
    small_drive = small.Scat*small.input_projection*Ls;
    outer_maps = @(ids) hybrid_maps(ids,Po,Sd,Rd,Ps, ...
        small.output_interpolation,small_drive,zk);
    outer = build_interface(Po,outer_maps,tol,leaf_points,'Outer sphere');
    [BoI,CoI] = outer_maps(outer.nodes);
    [outer.Scat,errors.outer_selected_residual] = ...
        selected_response(H,FH,BoI,CoI);
    outer.zk = zk;
    outer.alpha = alpha;
    timings.outer_interface = toc(stage);
    fprintf('Outer scattering matrix: %d by %d\n',size(outer.Scat));

    % Held-out analytic fields have no sources inside the outer sphere.
    % No BIE solutions or training waves were used for either ID.
    stage = tic;
    [incd,labels] = incident_fields(Sd.r,zk,geometry.outer_center);
    incs = incident_fields(Ss.r,zk,geometry.outer_center);
    din_s = Ls*incident_fields(Ps.r,zk,geometry.outer_center);
    din_o = Lo*incident_fields(Po.r,zk,geometry.outer_center);
    rhs_full = [Rd*incd;Rs*incs];
    rhs_hybrid = [Rd*incd;small.Scat*(small.input_projection*din_s)];
    zh = FH\rhs_hybrid;
    errors.hybrid_solve = relative_errors(H*zh,rhs_hybrid);
    zh_outer = FH\(BoI*(outer.input_projection*din_o));
    bout = outer.Scat*(outer.input_projection*din_o);
    dout = outer.output_interpolation*bout;
    errors.outer_incoming_rhs = relative_errors( ...
        BoI*(outer.input_projection*din_o),rhs_hybrid);

    % Full reference uses the actual NRCCIE cross blocks and near corrections.
    Ads = physical_matrix(Ss,Sd,opts,eps_quad);
    Asd = physical_matrix(Sd,Ss,opts,eps_quad);
    A_full = [Ad,Ads;Asd,As];
    F_full = decomposition(A_full,'lu');
    qref = F_full\rhs_full;
    errors.full_bie_residual = relative_errors(A_full*qref,rhs_full);
    qd_hybrid = zh(1:nd,:);
    small_incident_coefficients = small.input_projection*din_s+Gsd*qd_hybrid;
    % Recover the small physical density only for diagnostics, not construction.
    qs_hybrid = Fs\(BsI*small_incident_coefficients);
    q_hybrid = [qd_hybrid;qs_hybrid];
    errors.hybrid_full_bie_residual = relative_errors(A_full*q_hybrid,rhs_full);
    errors.hybrid_density = relative_errors(q_hybrid,qref);
    errors.small_recovery = relative_errors(CsI*qs_hybrid,zh(nd+1:end,:));

    % -R*K coupling omits the charge-continuity term of general NRCCIE
    % densities. Quantify it using solved physical fields rather than assuming
    % it vanishes on every column of a kernel map.
    Ksd = nrccie_matrix(Ss,Sd,zk);
    Kds = nrccie_matrix(Sd,Ss,zk);
    errors.small_cross_compatibility = relative_errors( ...
        -Rd*Ksd*qref(nd+1:end,:),Ads*qref(nd+1:end,:));
    errors.large_cross_compatibility = relative_errors( ...
        -Rs*Kds*qref(1:nd,:),Asd*qref(1:nd,:));
    A_field = [Ad,-Rd*Ksd;-Rs*Kds,As];
    q_field = A_field\rhs_full;
    errors.field_coupled_full_residual = relative_errors(A_full*q_field,rhs_full);
    clear F_full A_field

    errors.small_incoming_rhs = relative_errors( ...
        BsI*small.input_projection*din_s,Rs*incs);
    errors.large_to_small_rhs = relative_errors(BsI*Gsd*qref(1:nd,:), ...
        Rs*Kds*qref(1:nd,:));
    errors.small_to_large_rhs = relative_errors(Gds*CsI*qref(nd+1:end,:), ...
        Rd*Ksd*qref(nd+1:end,:));

    proxy_hybrid = complex(zeros(4*Po.npts,numel(labels)));
    proxy_ref = proxy_hybrid;
    proxy_outer_direct = proxy_hybrid;
    for first = 1:leaf_points:Po.npts
      ids = first:min(first+leaf_points-1,Po.npts);
      rows = reshape(4*(ids-1)+(1:4).',1,[]);
      [~,Co] = hybrid_maps(ids,Po,Sd,Rd,Ps, ...
          small.output_interpolation,small_drive,zk,false);
      proxy_hybrid(rows,:) = Co*zh;
      proxy_outer_direct(rows,:) = Co*zh_outer;
      T = sphere_subset(Po,ids);
      proxy_ref(rows,:) = maxwell_trace_matrix(T)*( ...
          nrccie_matrix(Sd,T,zk)*qref(1:nd,:) ...
          +nrccie_matrix(Ss,T,zk)*qref(nd+1:end,:));
    end
    errors.outer_output_id = relative_errors( ...
        outer.output_interpolation*(CoI*zh_outer),proxy_outer_direct);
    errors.outer_proxy_vs_hybrid = relative_errors(dout,proxy_hybrid);
    errors.outer_proxy_vs_full = relative_errors(dout,proxy_ref);
    clear Co proxy_outer_direct

    direction = directions(48,0.37);
    target_radii = [1.6,2.4,4.0];
    targets.r = geometry.outer_center+[target_radii(1)*direction, ...
        target_radii(2)*direction,target_radii(3)*direction];
    Qt = stratton_chu_matrix(Po,targets,zk);
    small_target = stratton_chu_matrix(Ps,targets,zk)*small.output_interpolation;
    Kdt = nrccie_matrix(Sd,targets,zk);
    Kst = nrccie_matrix(Ss,targets,zk);
    field_ref = Kdt*qref(1:nd,:)+Kst*qref(nd+1:end,:);
    field_hybrid = Kdt*zh(1:nd,:)+small_target*zh(nd+1:end,:);
    field_outer = (Qt*outer.output_interpolation)*bout;
    field_fieldcoupled = Kdt*q_field(1:nd,:)+Kst*q_field(nd+1:end,:);
    errors.hybrid_target_vs_full = relative_errors(field_hybrid,field_ref);
    errors.outer_target_vs_full = relative_errors(field_outer,field_ref);
    errors.outer_target_vs_hybrid = relative_errors(field_outer,field_hybrid);
    errors.outer_green_hybrid = relative_errors(Qt*proxy_hybrid,field_hybrid);
    errors.outer_green_full = relative_errors(Qt*proxy_ref,field_ref);
    errors.field_coupled_target_vs_full = relative_errors(field_fieldcoupled,field_ref);
    errors.electric_target = component_errors(field_outer,field_ref,1:3);
    errors.magnetic_target = component_errors(field_outer,field_ref,4:6);
    for j = 1:numel(target_radii)
      rows = 6*(j-1)*size(direction,2)+(1:6*size(direction,2));
      errors.target_by_radius(j,:) = relative_errors(field_outer(rows,:),field_ref(rows,:));
    end
    % Dropping the feedback on the small sphere provides a control: agreement
    % must come from coupling, not from a negligibly interacting geometry.
    b_no_feedback = small.Scat*(small.input_projection*din_s);
    qd_no_feedback = Ad\(Rd*incd+Gds*b_no_feedback);
    field_no_feedback = Kdt*qd_no_feedback+small_target*b_no_feedback;
    errors.no_feedback_vs_full = relative_errors(field_no_feedback,field_ref);
    timings.validation = toc(stage);

    errors.small_operator_b = small.errors.b;
    errors.small_operator_c = small.errors.c;
    errors.outer_operator_b = outer.errors.b;
    errors.outer_operator_c = outer.errors.c;
    stats = struct('small_physical_dofs',ns,'large_physical_dofs',nd, ...
        'full_system_size',nd+ns,'hybrid_system_size',nd+ms, ...
        'small_sphere_points',Ps.npts,'small_skeleton_points',numel(small.nodes), ...
        'outer_sphere_points',Po.npts,'outer_skeleton_points',numel(outer.nodes), ...
        'small_scattering_size',size(small.Scat),'outer_scattering_size',size(outer.Scat));
    timings.total = toc(total_timer);
    accuracy_limit = max(1e-4,100*tol);
    if norder<12, accuracy_limit = max(accuracy_limit,5e-3); end
    passed = all(isfinite(field_outer(:))) ...
        && errors.small_selected_residual<1e-10 ...
        && errors.outer_selected_residual<1e-10 ...
        && max(errors.full_bie_residual)<1e-10 ...
        && max(errors.hybrid_solve)<1e-10 ...
        && max([small.errors.b,small.errors.c,outer.errors.b,outer.errors.c])<10*tol ...
        && numel(small.nodes)<Ps.npts && numel(outer.nodes)<Po.npts ...
        && max(errors.outer_target_vs_full)<accuracy_limit ...
        && max(errors.hybrid_full_bie_residual)<accuracy_limit;
    result = struct('passed',passed,'small',small,'outer',outer, ...
        'geometry',geometry,'norder',norder,'nproxy',nproxy,'tol',tol, ...
        'zk',zk,'alpha',alpha,'eps_quad',eps_quad,'stats',stats, ...
        'errors',errors,'timings',timings,'labels',{labels}, ...
        'targets',targets,'target_radii',target_radii,'field_ref',field_ref, ...
        'field_hybrid',field_hybrid,'field_outer',field_outer, ...
        'accuracy_limit',accuracy_limit,'physical_reference', ...
        'Same-discretization full two-body corrected NRCCIE; not a mesh-convergence claim');
    out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if ~isfolder(out_dir), mkdir(out_dir); end
    tag = sprintf('nested_ellipsoid_p%d_q%d_%s',norder,nproxy, ...
        strrep(sprintf('%.0e',tol),'-','m'));
    result.output_file = fullfile(out_dir,[tag,'.mat']);
    result.figure_file = fullfile(out_dir,[tag,'.png']);
    save(result.output_file,'-struct','result','-v7.3');
    print_errors(errors);
    fprintf('Full/hybrid dimensions: %d / %d; elapsed %.1f s\n', ...
        nd+ns,nd+ms,timings.total);
    fprintf('Saved %s\n',result.output_file);
    plot_result(result);
    assert(passed,'Nested ellipsoid validation failed; inspect the saved errors.');
    fprintf('NESTED_ELLIPSOID_SCATTERING_PASS\n');
end

function A = physical_matrix(S,T,opts,eps_quad)
    if isempty(T), nt = S.npts; else, nt = T.npts; end
    corrections = nrccie_quad_corr_block(S,eps_quad,opts.zk,T,128);
    A = nrccie_system_block(S,T,corrections,opts,1:3*nt,1:3*S.npts);
end

function [B,C] = body_maps(ids,P,S,R,zk)
    Q = sphere_subset(P,ids);
    B = -R*stratton_chu_matrix(Q,S,zk);
    C = maxwell_trace_matrix(Q)*nrccie_matrix(S,Q,zk);
end

function [B,C] = hybrid_maps(ids,P,S,R,Ps,Us,small_drive,zk,build_b)
    if nargin<9, build_b = true; end
    Q = sphere_subset(P,ids);
    if build_b
      B = [-R*stratton_chu_matrix(Q,S,zk); ...
           -small_drive*stratton_chu_matrix(Q,Ps,zk)];
    else
      B = [];
    end
    L = maxwell_trace_matrix(Q);
    C = [L*nrccie_matrix(S,Q,zk), L*stratton_chu_matrix(Ps,Q,zk)*Us];
end

function P = sphere_subset(P,ids)
    for name = {'r','n','dru','drv'}
      P.(name{1}) = P.(name{1})(:,ids);
    end
    P.wts = P.wts(ids);
    P.npts = numel(ids);
end

function interface = build_interface(P,mapfun,tol,leaf_points,label)
    % Stream the normalization and the final operator checks, retaining only
    % a single batch of B columns / C rows at a time.
    bnorm2 = 0;
    cnorm2 = 0;
    for first = 1:leaf_points:P.npts
      [B,C] = mapfun(first:min(first+leaf_points-1,P.npts));
      bnorm2 = bnorm2+norm(B,'fro')^2;
      cnorm2 = cnorm2+norm(C,'fro')^2;
    end
    bnorm = sqrt(bnorm2);
    cnorm = sqrt(cnorm2);
    assert(bnorm>0 && cnorm>0);
    fprintf('%s: streamed norms; starting multilevel joint ID\n',label);
    rowfun = @(ids) joint_rows(ids,mapfun,bnorm,cnorm);
    [nodes,U,info] = maxwell_multilevel_row_id(P,rowfun,tol,leaf_points);
    [BI,CI] = mapfun(nodes);
    berror2 = 0;
    cerror2 = 0;
    for first = 1:leaf_points:P.npts
      ids = first:min(first+leaf_points-1,P.npts);
      dofs = reshape(4*(ids-1)+(1:4).',1,[]);
      [B,C] = mapfun(ids);
      berror2 = berror2+norm(B-BI*U(dofs,:)','fro')^2;
      cerror2 = cerror2+norm(C-U(dofs,:)*CI,'fro')^2;
    end
    errors = struct('b',sqrt(berror2/bnorm2),'c',sqrt(cerror2/cnorm2));
    dofs = reshape(4*(nodes-1)+(1:4).',1,[]);
    assert(norm(U(dofs,:)-eye(numel(dofs)),'fro')<1e-11);
    interface = struct('proxy',P,'nodes',nodes,'dofs',dofs, ...
        'skeleton',sphere_subset(P,nodes),'input_projection',U', ...
        'output_interpolation',U,'id_info',info,'errors',errors, ...
        'normalization',struct('b_frobenius',bnorm,'c_frobenius',cnorm), ...
        'trace_order','Eu,Ev,Hu,Hv per node; unit tangents; no weights');
    fprintf('%s: %d/%d points; B/C errors %.3e / %.3e\n', ...
        label,numel(nodes),P.npts,errors.b,errors.c);
end

function F = joint_rows(ids,mapfun,bnorm,cnorm)
    [B,C] = mapfun(ids);
    F = [B'/bnorm,C/cnorm];
end

function [Scat,residual] = selected_response(A,F,B,C)
    Scat = complex(zeros(size(C,1),size(B,2)));
    residual2 = 0;
    for first = 1:128:size(B,2)
      ids = first:min(first+127,size(B,2));
      X = F\B(:,ids);
      Scat(:,ids) = C*X;
      residual2 = residual2+norm(A*X-B(:,ids),'fro')^2;
    end
    residual = sqrt(residual2)/norm(B,'fro');
end

function [EH,labels] = incident_fields(r,zk,center)
    d = directions(6,0.73);
    EH = complex(zeros(6*size(r,2),16));
    labels = cell(1,16);
    for j = 1:6
      [~,axis_id] = min(abs(d(:,j)));
      axis_vector = zeros(3,1); axis_vector(axis_id) = 1;
      u = cross(d(:,j),axis_vector); u = u/norm(u);
      v = cross(d(:,j),u);
      phase = exp(1i*zk*(d(:,j).'*r));
      EH(:,2*j-1) = reshape([u*phase;v*phase],[],1);
      EH(:,2*j) = reshape([v*phase;-u*phase],[],1);
      labels{2*j-1} = sprintf('plane_%d_u',j);
      labels{2*j} = sprintf('plane_%d_v',j);
    end
    for j = 1:2
      EH(:,12+j) = (EH(:,2*j-1)+1i*EH(:,2*j))/sqrt(2);
      labels{12+j} = sprintf('circular_%d',j);
    end
    source.r = center+[0.2;-0.3;3];
    source.edips = [1;0.3i;0.2];
    source.hdips = [0.2;1;0.1i];
    [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'ed');
    EH(:,15) = reshape([E;H],[],1); labels{15} = 'electric_dipole';
    [E,H] = em3d.incoming_sources(zk,source,struct('r',r),'hd');
    EH(:,16) = reshape([E;H],[],1); labels{16} = 'magnetic_dipole';
end

function d = directions(n,angle)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5))+angle;
    s = sqrt(1-z.^2);
    d = [s.*cos(phi);s.*sin(phi);z];
end

function e = relative_errors(value,reference)
    e = vecnorm(value-reference,2,1)./max(vecnorm(reference,2,1),realmin);
end

function e = component_errors(value,reference,components)
    rows = reshape(6*(0:size(value,1)/6-1)+components(:),1,[]);
    e = relative_errors(value(rows,:),reference(rows,:));
end

function print_errors(errors)
    names = fieldnames(errors);
    fprintf('\nMaximum relative errors over independent incident fields:\n');
    for j = 1:numel(names)
      value = errors.(names{j});
      fprintf('  %-36s %.6e\n',names{j},max(value(:)));
    end
end

function plot_result(result)
    fig = figure('Visible','off','Color','w','Position',[80,80,1250,520]);
    layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(layout);
    hold(ax,'on');
    t = linspace(0,2*pi,300);
    g = result.geometry;
    plot(ax,g.outer_center(1)+g.outer_radius*cos(t), ...
        g.outer_center(2)+g.outer_radius*sin(t),'k-','LineWidth',1.5);
    plot(ax,g.small_center(1)+g.small_radius*cos(t), ...
        g.small_center(2)+g.small_radius*sin(t),'Color',[0.1,0.6,0.25],'LineWidth',1.5);
    fill(ax,g.large_center(1)+g.large_axes(1)*cos(t), ...
        g.large_center(2)+g.large_axes(2)*sin(t),[0.4,0.6,0.85]);
    fill(ax,g.small_center(1)+g.small_axes(1)*cos(t), ...
        g.small_center(2)+g.small_axes(2)*sin(t),[0.95,0.6,0.2]);
    axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'x / wavelength'); ylabel(ax,'y / wavelength');
    title(ax,'Equatorial section: nested interfaces, disjoint bodies');
    legend(ax,{'Outer sphere','Small sphere','Large physical BIE','Small compressed body'}, ...
        'Location','southoutside');
    ax = nexttile(layout);
    semilogy(ax,1:numel(result.labels),result.errors.outer_target_vs_full,'o-', ...
        1:numel(result.labels),result.errors.hybrid_target_vs_full,'s-', ...
        1:numel(result.labels),result.errors.no_feedback_vs_full,'x-','LineWidth',1.2);
    grid(ax,'on'); xlabel(ax,'Independent incident field'); ylabel(ax,'Relative E/H error');
    title(ax,sprintf('Full NRCCIE comparison | small %d, outer %d DOFs', ...
        size(result.small.Scat,1),size(result.outer.Scat,1)));
    legend(ax,{'Outer compressed scattering','Hybrid coupled system', ...
        'Feedback on small body omitted'},'Location','southoutside');
    exportgraphics(fig,result.figure_file,'Resolution',150);
    close(fig);
end
