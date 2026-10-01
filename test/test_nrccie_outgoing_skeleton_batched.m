% Outgoing field skeleton on the 7,986-node analytic-test ellipsoid.
% Connected source blocks share one closed proxy; no full-source K or P.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(29)

semi_axes = [1.2,0.85,0.65];
patches_per_axis = [3,3,4];
surface_order = 10;
zk = 0.9;
tol = 1e-9;
max_patches = 6;
proxy_radius = 4;
validation_radius = 2*proxy_radius;
proxy_order = 6;
nvalidation = 128;
output_dir = fullfile(root,'diagnostics','outgoing_skeleton_batched');
if ~isfolder(output_dir), mkdir(output_dir); end
diary(fullfile(output_dir,'run.log'))
S = geometries.ellipsoid(semi_axes,patches_per_axis,[0;0;0],surface_order,11);
T = geometries.ellipsoid(proxy_radius*[1,1,1],[1,1,1],[0;0;0],proxy_order,11);
assert(S.npts==7986 && min(vecnorm(T.r))>max(vecnorm(S.r)))
fprintf('Source: %d patches, %d nodes; closed proxy: radius %.2f, %d nodes\n', ...
    S.npatches,S.npts,proxy_radius,T.npts)
construction_timer = tic;
[Jout,blocks,Pmerge,info] = hierarchical_nrccie_outgoing_skeleton( ...
    S,T,zk,tol,max_patches);
construction_seconds = toc(construction_timer);
skeleton_nodes = (Jout(1:3:end)+2)/3;
skeleton_r = S.r(:,skeleton_nodes);
nparts = numel(blocks);
assert(all(sum(info.adjacency,2)==4))
assert(all(info.part_patch_counts<=max_patches))
assert(isequal(sort([blocks.nodes]),1:S.npts))
assert(isequal(reshape(Jout,3,[]),3*skeleton_nodes+(-2:0).'))
assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
final_counts = accumarray(info.part_id(S.patch_id(skeleton_nodes)),1,[nparts,1]);
fprintf('%d edge-connected parts; construction %.2f s\n',nparts,construction_seconds)
fprintf('Part   Original nodes   Local skeleton   Final skeleton\n')
for part = 1:nparts
    fprintf('%4d   %14d   %14d   %14d\n',part,info.part_node_counts(part), ...
        info.local_skeleton_counts(part),final_counts(part))
end
fprintf('Merge: %d candidates -> %d complete nodes (%d DOFs)\n', ...
    numel(info.candidate_nodes),numel(skeleton_nodes),numel(Jout))

% Only the final skeleton is evaluated all at once.
Sskel = node_subset(S,skeleton_nodes);
Kskel = nrccie_matrix(Sskel,T,zk);
proxy_w = sqrt(repelem(T.wts(:),6));
rotations = sample_so3(nvalidation);
directions = reshape(rotations(:,3,:),3,[]);
validation_target.r = validation_radius*directions;
Kvalidation_skel = nrccie_matrix(Sskel,validation_target,zk);
identity = eye(numel(Jout));
proxy_squared = 0;
proxy_norm_squared = 0;
validation_squared = 0;
validation_norm_squared = 0;
source_w = sqrt(repelem(S.wts(:),3));
q = randn(3*S.npts,4)+1i*randn(3*S.npts,4);
q = q./vecnorm(q);
density = q./source_w;
equivalent_density = complex(zeros(numel(Jout),size(density,2)));
exact_field = complex(zeros(6*nvalidation,size(density,2)));
for part = 1:nparts
    nodes = blocks(part).nodes;
    source = node_subset(S,nodes);
    dofs = reshape(3*nodes+(-2:0).',1,[]);
    merge_part = Pmerge(:,blocks(part).candidate_dofs);
    Ppart = merge_part*blocks(part).P;
    [present,local_columns] = ismember(Jout,dofs);
    assert(norm(Ppart(:,local_columns(present))-identity(:,present),'fro')<1e-11)
    w = sqrt(repelem(source.wts(:),3));
    Kpart = nrccie_matrix(source,T,zk);
    residual = (proxy_w.*(Kpart-Kskel*Ppart))./w.';
    proxy_squared = proxy_squared+norm(residual,'fro')^2;
    proxy_norm_squared = proxy_norm_squared+norm((proxy_w.*Kpart)./w.','fro')^2;
    clear Kpart residual
    Kvalidation_part = nrccie_matrix(source,validation_target,zk);
    validation_squared = validation_squared+ ...
        norm((Kvalidation_part-Kvalidation_skel*Ppart)./w.','fro')^2;
    validation_norm_squared = validation_norm_squared+ ...
        norm(Kvalidation_part./w.','fro')^2;

    % Apply the two stored factors, independently of the explicit block fit.
    equivalent_density = equivalent_density+ ...
        merge_part*(blocks(part).P*density(dofs,:));
    exact_field = exact_field+Kvalidation_part*density(dofs,:);
    clear Kvalidation_part Ppart
end
proxy_error = sqrt(proxy_squared/proxy_norm_squared);
validation_error = sqrt(validation_squared/validation_norm_squared);
compressed_field = Kvalidation_skel*equivalent_density;
density_error = norm(compressed_field-exact_field,'fro')/norm(exact_field,'fro');
fprintf('Errors: proxy %.3e (bound %.3e), independent field at radius %.2f %.3e, density %.3e\n', ...
    proxy_error,info.proxy_error_bound,validation_radius,validation_error,density_error)
fprintf('Largest single proxy matrix %.2f MiB; full matrix would be %.2f MiB\n', ...
    info.largest_proxy_mib,info.full_proxy_mib)
fprintf('Stored interpolation factors: %.2f MiB\n',info.interpolation_mib)
assert(proxy_error<1.1*tol && proxy_error<=1.1*info.proxy_error_bound+200*eps)
assert(validation_error<3*tol && density_error<10*tol)
save(fullfile(output_dir,'skeleton_results.mat'), ...
    'semi_axes','patches_per_axis','surface_order','zk','tol','max_patches', ...
    'proxy_radius','proxy_order','validation_radius','nvalidation','Jout','skeleton_nodes','skeleton_r', ...
    'info','final_counts','construction_seconds','proxy_error','validation_error','density_error')

% Show connected parts from opposite sides, then the proxy and final points.
views = [35,24;215,-24];
for kind = 1:2
    fig = figure('Color','w','Position',[100,100,1200,540]);
    layout = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    for panel = 1:2
        ax = nexttile(layout);
        if kind==1
            h = plot(S,info.part_id);
            colors = turbo(nparts);
            colormap(ax,colors([1:2:end,2:2:end],:))
            clim(ax,[0.5,nparts+0.5])
            set(h,'EdgeColor','none','FaceAlpha',1)
            hold(ax,'on')
            wireframe(S,struct('nfac',0.2))
            view(ax,views(panel,:))
        else
            if panel==1
                proxy_handle = plot(T,zeros(T.npatches,1));
                set(proxy_handle,'FaceAlpha',0.08,'EdgeColor','none')
                hold(ax,'on')
                scatter3(ax,T.r(1,:),T.r(2,:),T.r(3,:),5,[0.15,0.5,0.9],'filled')
            end
            h = plot(S,zeros(S.npatches,1));
            set(h,'FaceColor',[0.86,0.89,0.93],'Vertices',0.985*h.Vertices, ...
                'FaceAlpha',1,'EdgeColor','none')
            hold(ax,'on')
            if panel==1
                set(proxy_handle,'FaceColor',[0.15,0.5,0.9])
            else
                wireframe(S,struct('nfac',0.2))
            end
            scatter3(ax,skeleton_r(1,:),skeleton_r(2,:),skeleton_r(3,:), ...
                24,[0.85,0.12,0.10],'filled','MarkerEdgeColor','w','LineWidth',0.3)
            view(ax,views(1,:))
        end
        axis(ax,'equal')
        axis(ax,'padded')
        axis(ax,'off')
        camproj(ax,'orthographic')
        hold(ax,'off')
    end
    if kind==1, plot_name = 'connected_parts'; else, plot_name = 'proxy_and_skeleton'; end
    exportgraphics(fig,fullfile(output_dir,[plot_name,'.png']),'Resolution',180)
    savefig(fig,fullfile(output_dir,[plot_name,'.fig']))
end
disp('NRCCIE_OUTGOING_SKELETON_BATCHED_PASS')
diary off

function source = node_subset(S,nodes)
source.npts = numel(nodes);
source.r = S.r(:,nodes);
source.dru = S.dru(:,nodes);
source.drv = S.drv(:,nodes);
source.wts = S.wts(nodes);
end
