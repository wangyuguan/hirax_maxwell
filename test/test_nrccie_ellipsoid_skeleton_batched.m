% Compress the original analytic-test ellipsoid in connected patch groups.
% Never assemble the full-surface training RHS or interpolation matrix.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(23)

semi_axes = [1.2,0.85,0.65];
patches_per_axis = [3,3,4];
surface_order = 10;
zk = 0.9;
alpha = 1;
tol = 1e-9;
nsample = 512;
nvalidation = 64;
max_patches = 6;
output_dir = fullfile(root,'diagnostics','ellipsoid_skeleton_batched');
if ~isfolder(output_dir), mkdir(output_dir); end
diary(fullfile(output_dir,'run.log'))

S = geometries.ellipsoid( ...
    semi_axes,patches_per_axis,[0;0;0],surface_order,11);
assert(S.npts==7986)
rotations = sample_so3(nsample);
validation_rotations = sample_so3(nvalidation);
construction_timer = tic;
[Jin,blocks,Pmerge,info] = hierarchical_nrccie_incoming_skeleton( ...
    S,zk,alpha,rotations,tol,max_patches);
part_id = info.part_id;
adjacency = info.adjacency;
nparts = numel(blocks);
assert(all(sum(adjacency,2)==4), ...
    'Every quadrilateral on this closed ellipsoid must have four neighbors.')
assert(all(conncomp(graph(adjacency))==1))
part_patch_counts = info.part_patch_counts;
assert(all(part_patch_counts<=max_patches))
assert(isequal(sort([blocks.nodes]),1:S.npts))
fprintf('Ellipsoid: %d patches, %d nodes, %d DOFs; k %.2f\n', ...
    S.npatches,S.npts,3*S.npts,zk)
fprintf('%d connected parts, at most %d patches per part; connectivity PASS\n', ...
    nparts,max_patches)

candidate_nodes = info.candidate_nodes;
part_node_counts = info.part_node_counts;
local_skeleton_counts = info.local_skeleton_counts;
local_errors = info.local_errors;
for part = 1:nparts
    fprintf('Part %2d/%d: %d patches, %4d nodes -> %3d candidates, error %.2e\n', ...
        part,nparts,part_patch_counts(part),part_node_counts(part), ...
        local_skeleton_counts(part),local_errors(part))
end
assert(numel(unique(candidate_nodes))==numel(candidate_nodes))
merge_error = info.merge_error;
skeleton_nodes = (Jin(1:3:end)+2)/3;
skeleton_r = S.r(:,skeleton_nodes);
assert(isequal(reshape(Jin,3,[]),3*skeleton_nodes+(-2:0).'))
assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
fprintf('Merge: %d candidate nodes -> %d final nodes, error %.2e\n', ...
    numel(candidate_nodes),numel(skeleton_nodes),merge_error)
% Apply the returned local/merge factors to training and independent waves.
T = node_subset(S,skeleton_nodes);
Bskel = nrccie_planewave_rhs(T,zk,alpha,rotations);
Bvalidation_skel = nrccie_planewave_rhs(T,zk,alpha,validation_rotations);
candidate_rhs = Pmerge*Bskel;
candidate_validation_rhs = Pmerge*Bvalidation_skel;
identity = eye(numel(Jin));
training_squared = 0;
training_norm_squared = 0;
validation_squared = zeros(1,2*nvalidation);
validation_norm_squared = zeros(1,2*nvalidation);
for part = 1:nparts
    nodes = blocks(part).nodes;
    T = node_subset(S,nodes);
    Bpart = nrccie_planewave_rhs(T,zk,alpha,rotations);
    candidate_dofs = blocks(part).candidate_dofs;
    Ppart = blocks(part).P*Pmerge(candidate_dofs,:);
    global_dofs = reshape(3*nodes+(-2:0).',1,[]);
    [present,local_rows] = ismember(Jin,global_dofs);
    assert(norm(Ppart(local_rows(present),:)-identity(present,:),'fro')<1e-11)
    w = sqrt(repelem(T.wts(:),3));
    reconstructed = blocks(part).P*candidate_rhs(candidate_dofs,:);
    training_squared = training_squared+norm(w.*(Bpart-reconstructed),'fro')^2;
    training_norm_squared = training_norm_squared+norm(w.*Bpart,'fro')^2;
    clear Bpart
    Bvalidation = nrccie_planewave_rhs(T,zk,alpha,validation_rotations);
    reconstructed = blocks(part).P*candidate_validation_rhs(candidate_dofs,:);
    residual = w.*(Bvalidation-reconstructed);
    validation_squared = validation_squared+sum(abs(residual).^2,1);
    validation_norm_squared = validation_norm_squared+sum(abs(w.*Bvalidation).^2,1);
    clear Ppart Bvalidation residual
end
skeleton_seconds = toc(construction_timer);
training_error = sqrt(training_squared/training_norm_squared);
validation_error = sqrt(sum(validation_squared)/sum(validation_norm_squared));
worst_wave_error = max(sqrt(validation_squared./validation_norm_squared));
full_rhs_mib = info.full_rhs_mib;
largest_rhs_mib = info.largest_rhs_mib;
fprintf('Retained nodes: %d / %d (%.2f%%); retained DOFs: %d\n', ...
    numel(skeleton_nodes),S.npts,100*numel(skeleton_nodes)/S.npts,numel(Jin))
fprintf('Weighted errors: training %.3e, held-out %.3e, worst wave %.3e\n', ...
    training_error,validation_error,worst_wave_error)
fprintf('Largest single RHS: %.2f MiB; full RHS would be %.2f MiB\n', ...
    largest_rhs_mib,full_rhs_mib)
fprintf('Compression and streamed validation: %.2f s\n',skeleton_seconds)
fprintf('Composed training error bound: %.3e; stored factors: %.2f MiB\n', ...
    info.training_error_bound,info.interpolation_mib)
assert(training_error<tol, 'Global training error exceeds the requested tolerance.')
assert(training_error<=1.1*info.training_error_bound+200*eps)
assert(validation_error<3*tol && worst_wave_error<10*tol)
save(fullfile(output_dir,'skeleton_results.mat'), ...
    'semi_axes','patches_per_axis','surface_order','zk','alpha','tol', ...
    'nsample','nvalidation','max_patches','part_id','adjacency', ...
    'part_patch_counts','part_node_counts','local_skeleton_counts', ...
    'candidate_nodes','skeleton_nodes','skeleton_r','Jin','local_errors', ...
    'merge_error','training_error','validation_error','worst_wave_error', ...
    'full_rhs_mib','largest_rhs_mib','skeleton_seconds','info')

% Opposite views show both sides. Patch edges are the actual curved edges.
views = [35,24;215,-24];
for plot_kind = 1:2
    fig = figure('Color','w','Position',[100,100,1200,540]);
    layout = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    for side = 1:2
        ax = nexttile(layout);
        if plot_kind==1
            h = plot(S,part_id);
            colors = turbo(nparts);
            colormap(ax,colors([1:2:end,2:2:end],:))
            clim(ax,[0.5,nparts+0.5])
        else
            h = plot(S,zeros(S.npatches,1));
            % Inset only the display surface; markers keep their true coordinates.
            set(h,'FaceColor',[0.86,0.89,0.93],'Vertices',0.985*h.Vertices)
        end
        set(h,'EdgeColor','none','FaceAlpha',1)
        hold(ax,'on')
        wireframe(S,struct('nfac',0.2))
        if plot_kind==2
            scatter3(ax,skeleton_r(1,:),skeleton_r(2,:),skeleton_r(3,:), ...
                34,[0.85,0.12,0.10],'filled','MarkerEdgeColor','w','LineWidth',0.4)
        end
        axis(ax,'equal')
        axis(ax,'padded')
        axis(ax,'off')
        view(ax,views(side,:))
        camproj(ax,'orthographic')
        hold(ax,'off')
    end
    if plot_kind==1, plot_name = 'connected_parts'; else, plot_name = 'skeleton_points'; end
    exportgraphics(fig,fullfile(output_dir,[plot_name,'.png']),'Resolution',180)
    savefig(fig,fullfile(output_dir,[plot_name,'.fig']))
end
disp('NRCCIE_ELLIPSOID_SKELETON_BATCHED_PASS')
diary off

function T = node_subset(S,nodes)
T.npts = numel(nodes);
T.r = S.r(:,nodes);
T.n = S.n(:,nodes);
T.dru = S.dru(:,nodes);
T.drv = S.drv(:,nodes);
T.wts = S.wts(nodes);
end
