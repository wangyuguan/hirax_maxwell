%% Parameters
root = fileparts(mfilename('fullpath'));
geometry_file = fullfile(root,'data','common_local_geometry_order06.mat');
common_color = [0.64 0.79 0.91];
local_colors = [0.96 0.57 0.23; 0.29 0.69 0.48; ...
                0.68 0.49 0.82; 0.89 0.38 0.49];
show_mesh = true;
line_width = 0.3;
face_samples = 5;
edge_samples = 9;
local_xlim = [-25 25];
local_ylim = [-73 -49];
center_xlim = [-35 35];
center_ylim = [-35 35];
view_3d = [35 35];
z_scale = 1;  % Thickness multiplier for the 3D display only.

%% Load and sample the five meshes
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'));
clear pth dir
addpath(fullfile(root,'src'));
saved = load(geometry_file,'geometry','metadata');
geometry = saved.geometry;
common = sample_mesh(geometry.common,face_samples,edge_samples);
local = cell(1,4);
common_world = cell(1,4);
local_world = cell(1,4);
leaf_names = {'Upper left','Upper right','Lower left','Lower right'};
for k = 1:4
    local{k} = sample_mesh(geometry.local_tips{k},face_samples,edge_samples);
    R = geometry.rotation(:,:,k);
    shift = geometry.shift(:,k).';
    common_world{k} = common;
    common_world{k}.vertices = common.vertices*R.'+shift;
    common_world{k}.edges = common.edges*R.'+shift;
    local_world{k} = local{k};
    local_world{k}.vertices = local{k}.vertices*R.'+shift;
    local_world{k}.edges = local{k}.edges*R.'+shift;
end

%% Common part
fig_common = figure('Name','Common part','Color','w');
tiledlayout(fig_common,1,2,'TileSpacing','compact');
ax = nexttile;
draw_mesh(ax,common,common_color,show_mesh,line_width,1);
axis(ax,'equal'); axis(ax,'tight'); view(ax,2);
xlabel(ax,'u (mm)'); ylabel(ax,'v (mm)'); title(ax,'Common mesh');
ax = nexttile;
for k = 1:4
    draw_mesh(ax,common_world{k},common_color,show_mesh,line_width,1);
end
axis(ax,'equal'); axis(ax,'tight'); view(ax,2);
xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); title(ax,'Four copies');

%% Local parts
fig_local = figure('Name','Local parts','Color','w');
tiledlayout(fig_local,2,2,'TileSpacing','compact');
for k = 1:4
    ax = nexttile;
    draw_mesh(ax,local{k},local_colors(k,:),show_mesh,line_width,1);
    axis(ax,'equal'); view(ax,2);
    xlim(ax,local_xlim); ylim(ax,local_ylim);
    xlabel(ax,'u (mm)'); ylabel(ax,'v (mm)'); title(ax,leaf_names{k});
end

%% Final mesh
fig_final = figure('Name','Final mesh','Color','w','Position',[100 100 1300 800]);
layout = tiledlayout(fig_final,2,2,'TileSpacing','compact');
for panel = 1:3
    if panel == 1
        ax = nexttile(layout,1,[2 1]);
    else
        ax = nexttile(layout,2*(panel-1));
    end
    z_display = 1;
    if panel == 3, z_display = z_scale; end
    h = gobjects(1,5);
    for k = 1:4
        hc = draw_mesh(ax,common_world{k},common_color,show_mesh,line_width,z_display);
        if k == 1, h(1) = hc; end
        h(k+1) = draw_mesh(ax,local_world{k},local_colors(k,:),show_mesh,line_width,z_display);
    end
    axis(ax,'equal'); axis(ax,'tight'); view(ax,2);
    xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)');
    if panel == 1
        title(ax,'Final mesh');
        lg = legend(ax,h,[{'Common'},leaf_names],'Orientation','horizontal');
        lg.Layout.Tile = 'south';
    else
        xlim(ax,center_xlim); ylim(ax,center_ylim);
        title(ax,'Center');
    end
    if panel == 3
        view(ax,view_3d);
        ax.ZTick = [];
        title(ax,sprintf('3D: thickness %g mm, display x %g',saved.metadata.thickness,z_scale));
    end
end


function mesh = sample_mesh(S,nface,nedge)
% Sample faces for rendering and true patch perimeters for mesh lines.
[kinds,~,group] = unique([S.norders(:),S.iptype(:)],'rows');
ng = size(kinds,1);
F = cell(1,ng); B = cell(1,ng); E = cell(1,ng);
nv = zeros(1,ng); nf = nv; ne = nv;
for k = 1:ng
    p = kinds(k,1);
    if kinds(k,2) == 1
        uv = [0 1 0;0 0 1];
        F{k} = [1 2 3];
        B{k} = koorn.pols(p,uv);
        E{k} = koorn.pols(p,uv(:,[1 2 3 1]));
    else
        [u,v] = meshgrid(linspace(-1,1,nface));
        uv = [u(:).';v(:).'];
        index = reshape(1:nface^2,nface,nface);
        a = index(1:end-1,1:end-1); b = index(1:end-1,2:end);
        c = index(2:end,1:end-1); d = index(2:end,2:end);
        F{k} = [a(:),b(:),c(:);b(:),d(:),c(:)];
        s = linspace(-1,1,nedge); one = ones(size(s));
        perimeter = [s,one,-s,-one;-one,s,one,-s];
        B{k} = polytens.lege.pols(p,uv);
        E{k} = polytens.lege.pols(p,perimeter);
    end
    nv(k) = size(B{k},2); nf(k) = size(F{k},1); ne(k) = size(E{k},2)+1;
end
mesh.vertices = zeros(sum(nv(group)),3);
mesh.faces = zeros(sum(nf(group)),3);
mesh.edges = nan(sum(ne(group)),3);
vo = 0; fo = 0; eo = 0;
for j = 1:S.npatches
    k = group(j);
    coef = S.srccoefs{j}(1:3,:);
    mesh.vertices(vo+(1:nv(k)),:) = (coef*B{k}).';
    mesh.faces(fo+(1:nf(k)),:) = vo+F{k};
    mesh.edges(eo+(1:ne(k)-1),:) = (coef*E{k}).';
    vo = vo+nv(k); fo = fo+nf(k); eo = eo+ne(k);
end
end


function h = draw_mesh(ax,mesh,color,show_mesh,line_width,z_scale)
vertices = mesh.vertices;
vertices(:,3) = z_scale*vertices(:,3);
hold(ax,'on');
h = patch(ax,'Faces',mesh.faces,'Vertices',vertices, ...
    'FaceColor',color,'EdgeColor','none');
if show_mesh
    edges = mesh.edges;
    edges(:,3) = z_scale*edges(:,3);
    plot3(ax,edges(:,1),edges(:,2),edges(:,3), ...
        'Color',[0.18 0.22 0.27],'LineWidth',line_width,'HandleVisibility','off');
end
end
