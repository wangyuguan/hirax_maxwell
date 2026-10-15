function result = test_ellipsoid_shared_skeleton(tol,filename)
% TEST_ELLIPSOID_SHARED_SKELETON Use the same sphere nodes for input and output.
%
%  Input arguments:
%    * tol: joint block-row ID tolerance (1e-8)
%    * filename: full sphere scattering result (ellipsoid_p16_q20.mat)
%
%  Output arguments:
%    * result: common nodes, interpolation, square compressed Scat, and errors
%
%  Run startup.m first. This test reuses the saved full scattering matrix;
%  it does not assemble or solve another physical BIE. Every selected node
%  retains all four field values [Eu,Ev,Hu,Hv].

    if nargin < 1 || isempty(tol), tol = 1e-8; end
    output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
    if nargin < 2 || isempty(filename)
      filename = fullfile(output_dir,'ellipsoid_p16_q20.mat');
    end
    full = load(filename);
    P = full.proxy;
    zk = full.zk;
    L = maxwell_trace_matrix(P);
    nplane = 384;
    nsource = 96;
    fprintf('Training: %d plane waves and %d exterior dipole fields\n', ...
        2*nplane,6*nsource);

    % Incoming traces and their scattered responses span the two trace spaces.
    timer = tic;
    incoming = plane_traces(P,L,zk,directions(nplane,0));
    incoming = [incoming,dipole_traces(P,L,zk,2.5*directions(nsource,0.19))];
    outgoing = full.Scat*incoming;
    weights = repelem(P.wts(:),4);
    incoming = incoming./sqrt(sum(weights.*abs(incoming).^2,1));
    outgoing = outgoing./sqrt(sum(weights.*abs(outgoing).^2,1));
    [nodes,T,id_info] = maxwell_shared_row_id([incoming,outgoing],tol);
    dofs = reshape(4*nodes(:).'+(-3:0).',[],1);
    fprintf('Selected %d/%d common points (%d/%d scalar coordinates)\n', ...
        numel(nodes),P.npts,numel(dofs),4*P.npts);
    fprintf('Joint ID residual %.3e; largest interpolation entry %.3e\n', ...
        id_info.relative_residual,id_info.maxabsT);
    clear incoming outgoing

    % The restriction uses the same DOFs on both sides. T includes the lift;
    % taking full.Scat(dofs,dofs) alone would omit it.
    Scat = full.Scat(dofs,:)*T;
    build_time = toc(timer);
    Q = stratton_chu_matrix(P,full.targets,zk);
    target_matrix = Q*T;

    % New directions, both transverse polarizations, complex polarizations,
    % and dipoles at new positions. Original four BIE checks are also held out.
    incoming = full.incident_proxy;
    labels = full.labels;
    held_directions = directions(64,0.73);
    plane = plane_traces(P,L,zk,held_directions);
    incoming = [incoming,plane];
    labels = [labels,repmat({'new_plane'},1,size(plane,2))];
    circular = (plane(:,1:2:end)+1i*plane(:,2:2:end))/sqrt(2);
    incoming = [incoming,circular];
    labels = [labels,repmat({'new_complex_plane'},1,size(circular,2))];
    source_positions = directions(12,0.61).*repmat([2.5,3,4],1,4);
    dipoles = dipole_traces(P,L,zk,source_positions);
    incoming = [incoming,dipoles];
    labels = [labels,repmat({'new_dipole'},1,size(dipoles,2))];

    reference_trace = full.Scat*incoming;
    reference_target = Q*reference_trace;
    incident_skeleton = incoming(dofs,:);
    scattered_skeleton = Scat*incident_skeleton;
    scattered_trace = T*scattered_skeleton;
    scattered_target = target_matrix*scattered_skeleton;

    errors.incoming_trace = relative_errors(T*incident_skeleton,incoming);
    errors.outgoing_trace = relative_errors(T*reference_trace(dofs,:),reference_trace);
    errors.skeleton = relative_errors(scattered_skeleton,reference_trace(dofs,:));
    errors.composed_trace = relative_errors(scattered_trace,reference_trace);
    errors.target_compression = relative_errors(scattered_target,reference_target);
    errors.target_direct_bie = relative_errors(scattered_target(:,1:4),full.direct_target);
    errors.baseline_target_direct_bie = full.errors.target;

    fprintf('Held-out fields: %d\n',size(incoming,2));
    fprintf('Maximum incoming trace reconstruction: %.3e\n',max(errors.incoming_trace));
    fprintf('Maximum outgoing trace reconstruction: %.3e\n',max(errors.outgoing_trace));
    fprintf('Maximum compressed scattering trace error: %.3e\n',max(errors.composed_trace));
    fprintf('Maximum added target error from compression: %.3e\n',max(errors.target_compression));
    fprintf('Original BIE checks:     uncompressed       compressed\n');
    for j = 1:4
      fprintf('%-20s   %.3e          %.3e\n',full.labels{j}, ...
          full.errors.target(j),errors.target_direct_bie(j));
    end
    assert(isequal(size(Scat),[numel(dofs),numel(dofs)]));
    assert(norm(T(dofs,:)-eye(numel(dofs)),'fro')<1e-12);
    assert(numel(nodes)<P.npts,'No reduction in sphere node count.');
    assert(all(isfinite(Scat(:))),'Nonfinite compressed scattering matrix.');
    assert(id_info.relative_residual<tol,'Joint ID did not meet its training tolerance.');
    assert(max(errors.incoming_trace)<100*tol,'Held-out incoming trace check failed.');
    assert(max(errors.outgoing_trace)<100*tol,'Held-out outgoing trace check failed.');
    assert(max(errors.composed_trace)<100*tol,'Held-out scattering trace check failed.');
    assert(max(errors.target_compression)<100*tol, ...
        'Held-out target compression error exceeds 100 times the ID tolerance.');

    skeleton = struct('r',P.r(:,nodes),'n',P.n(:,nodes), ...
        'dru',P.dru(:,nodes),'drv',P.drv(:,nodes),'npts',numel(nodes));
    result = struct('Scat',Scat,'nodes',nodes,'dofs',dofs,'interpolation',T, ...
        'skeleton',skeleton,'proxy',P,'zk',zk,'tol',tol, ...
        'id_info',id_info,'source_file',filename,'build_time',build_time, ...
        'training_plane_directions',nplane,'training_dipole_positions',nsource, ...
        'training_source_radius',2.5,'labels',{labels},'errors',errors, ...
        'targets',full.targets,'target_matrix',target_matrix, ...
        'incident_skeleton',incident_skeleton,'scattered_skeleton',scattered_skeleton, ...
        'scattered_target',scattered_target,'reference_target',reference_target);
    tag = strrep(sprintf('%.0e',tol),'-','m');
    output_file = fullfile(output_dir,['ellipsoid_shared_skeleton_',tag,'.mat']);
    save(output_file,'-struct','result','-v7.3');
    fprintf('Saved %s\n',output_file);
    plot_result(result,output_dir,tag);
end

function F = plane_traces(P,L,zk,d)
    F = complex(zeros(4*P.npts,2*size(d,2)));
    for j = 1:size(d,2)
      [~,axis] = min(abs(d(:,j)));
      unit = zeros(3,1); unit(axis) = 1;
      u = cross(d(:,j),unit); u = u/norm(u);
      v = cross(d(:,j),u);
      phase = exp(1i*zk*(d(:,j).'*P.r));
      F(:,2*j-1) = L*reshape([u*phase;v*phase],[],1);
      F(:,2*j) = L*reshape([v*phase;-u*phase],[],1);
    end
end

function F = dipole_traces(P,L,zk,positions)
    F = complex(zeros(4*P.npts,6*size(positions,2)));
    for j = 1:size(positions,2)
      source.r = positions(:,j);
      for c = 1:3
        dipole = zeros(3,1); dipole(c) = 1;
        source.edips = dipole;
        source.hdips = dipole;
        [E,H] = em3d.incoming_sources(zk,source,P,'ed');
        F(:,6*j-6+c) = L*reshape([E;H],[],1);
        [E,H] = em3d.incoming_sources(zk,source,P,'hd');
        F(:,6*j-3+c) = L*reshape([E;H],[],1);
      end
    end
end

function d = directions(n,angle)
    z = 1-2*((0:n-1)+0.5)/n;
    phi = (0:n-1)*pi*(3-sqrt(5))+angle;
    s = sqrt(1-z.^2);
    d = [s.*cos(phi);s.*sin(phi);z];
    R = [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
    d = R*d;
end

function error = relative_errors(value,reference)
    error = vecnorm(value-reference,2,1)./vecnorm(reference,2,1);
end

function plot_result(result,output_dir,tag)
    fig = figure('Visible','off','Color','w','Position',[50,50,1120,440]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile;
    P = result.proxy;
    scatter3(P.r(1,:),P.r(2,:),P.r(3,:),7,[0.78,0.8,0.82],'filled');
    hold on;
    r = result.skeleton.r;
    scatter3(r(1,:),r(2,:),r(3,:),18,[0.1,0.4,0.65],'filled');
    axis equal; view(35,23); grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Same %d / %d sphere nodes for input and output', ...
        result.skeleton.npts,P.npts));
    legend({'Full sphere nodes','Shared skeleton nodes'}, ...
        'Location','southoutside','Orientation','horizontal');
    nexttile;
    e = result.errors;
    semilogy(e.incoming_trace,'-','LineWidth',1); hold on;
    semilogy(e.outgoing_trace,'-','LineWidth',1);
    semilogy(e.target_compression,'-','LineWidth',1);
    bounds = [e.incoming_trace,e.outgoing_trace,e.target_compression];
    ylim([0.7*min(bounds),2*max(bounds)]);
    xlabel('Held-out incident field'); ylabel('Relative error'); grid on;
    title(sprintf('Compressed scattering matrix: %d x %d',size(result.Scat)));
    legend({'Incoming trace','Outgoing trace','Added target error'}, ...
        'Location','best');
    output = fullfile(output_dir,['ellipsoid_shared_skeleton_',tag,'.png']);
    exportgraphics(fig,output,'Resolution',170);
    close(fig);
    fprintf('Saved %s\n',output);
end
