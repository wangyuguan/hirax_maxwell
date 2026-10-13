function [S,fig] = preview_dish(filename,visible)
%PREVIEW_DISH View an exported .go3 dish surface and its patch boundaries.
%   [S,fig] = preview_dish() loads the local order-4 f/D=0.23 rimfree mesh.
%   [S,fig] = preview_dish(filename,'off') renders without showing a window.
%   Coordinates retain the units of the input mesh. Figures are saved in
%   src/dish/results. Set DISH_MESH_ROOT to override the data/dish directory.

setup_dish();
dish_root = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(dish_root));
if nargin < 1 || isempty(filename)
    mesh_root = getenv('DISH_MESH_ROOT');
    if isempty(mesh_root)
        mesh_root = fullfile(project_root,'data','dish');
    end
    stem = 'dish_fd023_m_rimfree';
    filename = fullfile(mesh_root,stem,[stem,'_order4.go3']);
end
if nargin < 2, visible = 'on'; end
assert(isfile(filename), ...
    'Missing mesh %s. Supply a .go3 path or see src/dish/README.md.',filename);
S = surfer.load_from_file(filename);
[~,stem] = fileparts(filename);
fprintf('Loaded %s\n%d patches, %d nodes\n',filename,S.npatches,S.npts);

fig = figure('Name',['HIRAX dish: ',stem], ...
    'NumberTitle','off','Color','w','WindowStyle','normal', ...
    'Position',[100,150,1250,620],'Visible',visible);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('%s | %d patches',stem,S.npatches),'Interpreter','none');
ax1 = nexttile(layout);
surface_plot = plot(S,zeros(S.npatches,1));
set(surface_plot,'FaceColor',[0.35,0.62,0.82],'EdgeColor','none');
axis(ax1,'equal'); axis(ax1,'tight'); view(ax1,38,30);
xlabel(ax1,'x'); ylabel(ax1,'y'); zlabel(ax1,'z');
title(ax1,'Surface geometry');
camlight(ax1,'headlight'); lighting(ax1,'gouraud');

ax2 = nexttile(layout);
wireframe(S,struct('nfac',0));
set(findobj(ax2,'Type','line'),'LineWidth',0.35,'Color',[0.15,0.28,0.4]);
axis(ax2,'equal'); axis(ax2,'tight'); view(ax2,38,30); grid(ax2,'on');
xlabel(ax2,'x'); ylabel(ax2,'y'); zlabel(ax2,'z');
title(ax2,'High-order patch boundaries');
rotate3d(fig,'on');
drawnow;

output_dir = fullfile(dish_root,'results');
if ~isfolder(output_dir), mkdir(output_dir); end
output_stem = fullfile(output_dir,[stem,'_preview']);
savefig(fig,[output_stem,'.fig']);
exportgraphics(fig,[output_stem,'.png'],'Resolution',140);
fprintf('Preview ready: %s.fig\n',output_stem);
end
