% Outgoing skeleton, selected using a closed finite proxy surface.
% Small test only: no self matrix, singular quadrature, or Maxwell solve.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(29)

semi_axes = [1.2,0.85,0.65];
zk = 0.9;
tol = 1e-9;
proxy_radius = 4;
validation_radius = 2*proxy_radius;
proxy_order = 6;
nvalidation = 128;
S = geometries.ellipsoid(semi_axes,[1,1,1],[0;0;0],5,11);
T = geometries.ellipsoid(proxy_radius*[1,1,1],[1,1,1],[0;0;0],proxy_order,11);
assert(min(vecnorm(T.r))>max(vecnorm(S.r)))
Kproxy = nrccie_matrix(S,T,zk);
source_w = sqrt(repelem(S.wts(:),3));
proxy_w = sqrt(repelem(T.wts(:),6));
[Jout,Pout,proxy_error] = nrccie_outgoing_skeleton(S,proxy_w.*Kproxy,tol);
skeleton_nodes = (Jout(1:3:end)+2)/3;
assert(isequal(size(Pout),[numel(Jout),3*S.npts]))
assert(isequal(reshape(Jout,3,[]),3*skeleton_nodes+(-2:0).'))
assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
assert(norm(Pout(:,Jout)-eye(numel(Jout)),'fro')<1e-12)
measured_error = norm((proxy_w.*(Kproxy-Kproxy(:,Jout)*Pout))./source_w.','fro') ...
    /norm((proxy_w.*Kproxy)./source_w.','fro');
assert(abs(measured_error-proxy_error)<1e-13 && proxy_error<1.1*tol)

% Independent exterior targets do not participate in skeleton selection.
rotations = sample_so3(nvalidation);
directions = reshape(rotations(:,3,:),3,[]);
validation_targets.r = validation_radius*directions;
Kvalidation = nrccie_matrix(S,validation_targets,zk);
validation_error = norm((Kvalidation-Kvalidation(:,Jout)*Pout)./source_w.','fro') ...
    /norm(Kvalidation./source_w.','fro');
q = randn(3*S.npts,4)+1i*randn(3*S.npts,4);
q = q./vecnorm(q);
density = q./source_w;
equivalent_density = Pout*density;
exact_field = Kvalidation*density;
compressed_field = Kvalidation(:,Jout)*equivalent_density;
density_error = norm(compressed_field-exact_field,'fro') ...
    /norm(exact_field,'fro');
fprintf('Source: %d nodes; closed proxy: radius %.2f, %d nodes\n', ...
    S.npts,proxy_radius,T.npts)
fprintf('Outgoing skeleton: %d complete nodes, %d DOFs\n', ...
    numel(skeleton_nodes),numel(Jout))
fprintf('Errors: proxy %.3e, independent exterior field %.3e, density application %.3e\n', ...
    proxy_error,validation_error,density_error)
assert(validation_error<3*tol && density_error<10*tol)

% Cross-check proxy E,H against the existing NRCCIE tangential trace.
fields = reshape(Kproxy*density(:,1),6,T.npts);
E = fields(1:3,:);
H = fields(4:6,:);
alpha = 1;
tangent = -cross(T.n,H,1)+alpha*cross(T.n,cross(T.n,E,1),1);
expected_trace = [sum(T.dru.*tangent,1);sum(T.drv.*tangent,1)];
opts = struct('zk',zk,'alpha',alpha,'fmm',false);
trace = nrccie_apply(S,reshape(density(:,1),3,S.npts),T,[],opts);
assert(norm(trace(1:2,:)-expected_trace,'fro')<1e-12*norm(expected_trace,'fro'))
clear Kproxy

% Exact low-rank case with unequal area weights and pivots in two nodes.
Ssmall.wts = [1;4;9];
scaled_kernel = [4,.2,0,.1,0,.1i,.15,.1,.25i; ...
    0,.1i,.1,.2,1,.1,.1,.15i,.1];
small_w = sqrt(repelem(Ssmall.wts,3)).';
for nrow = 1:2
    Ksmall = scaled_kernel(1:nrow,:).*small_w;
    [Jsmall,Psmall,small_error] = nrccie_outgoing_skeleton(Ssmall,Ksmall,1e-12);
    assert(isequal(Jsmall,1:3*nrow))
    assert(norm(Psmall(:,Jsmall)-eye(numel(Jsmall)),'fro')<1e-12)
    assert(small_error<1e-12)
    d = randn(9,1)+1i*randn(9,1);
    assert(norm(Ksmall*d-Ksmall(:,Jsmall)*(Psmall*d))<1e-12*norm(Ksmall*d))
end
[Jzero,Pzero,zero_error] = nrccie_outgoing_skeleton(Ssmall,zeros(2,9),1e-12);
assert(isempty(Jzero) && isequal(size(Pzero),[0,9]) && zero_error==0)

output_dir = fullfile(root,'diagnostics','outgoing_skeleton');
if ~isfolder(output_dir), mkdir(output_dir); end
save(fullfile(output_dir,'skeleton_results.mat'), ...
    'semi_axes','zk','tol','proxy_radius','proxy_order','nvalidation','validation_radius', ...
    'Jout','Pout','skeleton_nodes','proxy_error','validation_error','density_error')
fig = figure('Color','w','Position',[100,100,1200,540]);
layout = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
for panel = 1:2
    ax = nexttile(layout);
    if panel==1
        h = plot(T,zeros(T.npatches,1));
        set(h,'FaceColor',[0.15,0.5,0.9],'FaceAlpha',0.08,'EdgeColor','none')
        hold(ax,'on')
        scatter3(ax,T.r(1,:),T.r(2,:),T.r(3,:),5,[0.15,0.5,0.9],'filled')
    end
    h = plot(S,zeros(S.npatches,1));
    set(h,'FaceColor',[0.86,0.89,0.93],'Vertices',0.985*h.Vertices, ...
        'FaceAlpha',1,'EdgeColor','none')
    hold(ax,'on')
    if panel==2, wireframe(S,struct('nfac',0.2)); end
    r = S.r(:,skeleton_nodes);
    scatter3(ax,r(1,:),r(2,:),r(3,:),24,[0.85,0.12,0.10], ...
        'filled','MarkerEdgeColor','w','LineWidth',0.3)
    axis(ax,'equal')
    axis(ax,'padded')
    axis(ax,'off')
    view(ax,35,24)
    camproj(ax,'orthographic')
    hold(ax,'off')
end
exportgraphics(fig,fullfile(output_dir,'proxy_and_skeleton.png'),'Resolution',180)
savefig(fig,fullfile(output_dir,'proxy_and_skeleton.fig'))
disp('NRCCIE_OUTGOING_SKELETON_PASS')
