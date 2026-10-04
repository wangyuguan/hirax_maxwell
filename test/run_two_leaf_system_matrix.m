% Dense off-diagonal NRCCIE block K=A21: source leaf 1 -> target leaf 2.
% Includes near quadrature corrections; no self jump. Each node has [ju;jv;rho].
% Column ID selects source components; plot their corresponding source nodes.
% Run on a server: even order 1 produces an 8856-by-8856 complex block.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Parameters: original shape, coarse panelization and order-1 quadrature
settings.leaf_ids = [1,2]; % TL/TR/BL/BR = 1/2/3/4
settings.leaf_radius = 69.25; % mm
settings.thickness = 0.6925;
settings.surface_gap = 4;
settings.norder = 1; % do not use 0: constant position patches degenerate
settings.zk = 2*pi/(10*settings.leaf_radius);
settings.alpha = 1;
settings.eps_quad = 1e-8;
settings.target_batch_nodes = 32;
settings.geometry_options = struct('norder',settings.norder, ...
    'chunkie_order',20,'chunkie_n0',3,'chunkie_nchs',3, ...
    'chunkie_newton_iterations',30, ...
    'rim_width',0.028128271246*settings.leaf_radius, ...
    'cap_collar_width',0.028128271246*settings.leaf_radius, ...
    'cap_mesh_spacing_center',0.5*settings.leaf_radius, ...
    'cap_mesh_spacing_side',0.5*settings.leaf_radius, ...
    'outline_refinement',1,'wall_profile_refinement',1);
assert(numel(settings.leaf_ids)==2 && numel(unique(settings.leaf_ids))==2 ...
    && all(ismember(settings.leaf_ids,1:4)))
assert(settings.norder>=1 && settings.norder==fix(settings.norder))
rng(79)

%% One original leaf shape, placed twice in the dish
[base,parts] = leaf_plate_surfer(settings.thickness,settings.geometry_options);
settings.geometry_options = parts.options;
center = parts.design.leaf_center_raw;
center(2) = center(2)+(settings.surface_gap-2.83*sqrt(2))/sqrt(2);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
leaves = cell(1,2);
for leaf = 1:2
    angle = angles(settings.leaf_ids(leaf));
    R = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
    settings.rotation(:,:,leaf) = R;
    settings.shift(:,leaf) = R*[center;0];
    leaves{leaf} = affine_transf(base,R,settings.shift(:,leaf));
end


C21 = nrccie_quad_corr_block(leaves{1},settings.eps_quad,settings.zk, ...
    leaves{2},settings.target_batch_nodes);
opts = struct('zk',settings.zk,'alpha',settings.alpha);
K = nrccie_system_block(leaves{1},leaves{2},C21,opts, ...
    1:3*leaves{2}.npts,1:3*leaves{1}.npts);


%% ID of the cross-leaf block
[Sk,RD,T] = id(K,1e-4);
save('id_result.mat','Sk','RD','T')

%% Plot the nodes represented by the selected components

selected_nodes = ceil(Sk(:)/3);
selected_components = mod(Sk(:)-1,3)+1; % 1=ju, 2=jv, 3=rho
skeleton_nodes = unique(selected_nodes);
r = leaves{1}.r;
rsk = r(:,skeleton_nodes);
fprintf('ID selected %d components on %d / %d source nodes.\n', ...
    numel(Sk),numel(skeleton_nodes),leaves{1}.npts)

fig = figure('Color','w','Position',[100,100,900,700]);
scatter3(r(1,:),r(2,:),r(3,:),6,[0.78,0.80,0.83],'filled')
hold on
scatter3(rsk(1,:),rsk(2,:),rsk(3,:),32,[0.85,0.12,0.10], ...
    'filled','MarkerEdgeColor','w','LineWidth',0.3)
axis equal
axis off
view(35,50)
hold off

output_dir = fullfile(root,'data');
if ~isfolder(output_dir), mkdir(output_dir); end
exportgraphics(fig,fullfile(output_dir,'two_leaf_id_nodes.png'),'Resolution',180)
savefig(fig,fullfile(output_dir,'two_leaf_id_nodes.fig'))
