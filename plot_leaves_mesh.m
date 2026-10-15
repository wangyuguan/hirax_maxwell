% View one, two and four copies of the current closed leaf geometry.
run(fullfile(fileparts(mfilename('fullpath')),'startup.m'))

opts = struct();
opts.norder = 4;
opts.wall_profile_refinement = 3;
opts.surface_gap = 4;
thickness = 0.6925;
[~,parts] = generate_leaves_surfer(thickness,opts);
original = parts.base_surfer;
original_parts = parts.base_parts;
original_leaves = parts.leaves;

% Exaggerate z in the axes without changing the geometry.
z_scale = 30;

figure('Name','One original leaf','Color','w');
ax = axes; hold(ax,'on')
upper_ids = [original_parts.top_core,original_parts.top_collar, ...
    original_parts.upper_wall];
lower_ids = [original_parts.lower_wall,original_parts.bottom_collar, ...
    original_parts.bottom_core];
surface_colors = ones(original.npatches,1);
surface_colors(lower_ids) = 2;
plot(original_leaves{1},surface_colors,'EdgeColor','none');
colormap(ax,[0.12 0.35 0.75;0.9 0.32 0.12]);
clim(ax,[1 2]);
axis(ax,'tight'); daspect(ax,[1 1 1/z_scale]); axis(ax,'off'); view(ax,35,25)
drawnow

figure('Name','Two original leaves','Color','w');
ax = axes; hold(ax,'on')
upper_leaf = slicesurfer(original_leaves{1},upper_ids);
lower_leaf = slicesurfer(original_leaves{1},lower_ids);
blue_surface = plot(original_leaves{2},zeros(original.npatches,1), ...
    'EdgeColor','none');
lower_surface = plot(lower_leaf,zeros(lower_leaf.npatches,1), ...
    'EdgeColor','none');
upper_surface = plot(upper_leaf,zeros(upper_leaf.npatches,1), ...
    'EdgeColor','none');
set(blue_surface,'FaceColor',[0.12 0.35 0.75],'EdgeColor','none')
set(lower_surface,'FaceColor',[0.20 0.65 0.40],'EdgeColor','none')
set(upper_surface,'FaceColor',[0.9 0.32 0.12],'EdgeColor','none')
axis(ax,'tight'); daspect(ax,[1 1 1/z_scale]); axis(ax,'off'); view(ax,35,25)
drawnow

figure('Name','Four original leaves','Color','w');
ax = axes; hold(ax,'on')
colors = [0.12 0.35 0.75;0.9 0.32 0.12; ...
    0.20 0.65 0.40;0.58 0.35 0.72];
surfaces = gobjects(1,4);
for leaf = 1:4
    surfaces(leaf) = plot(original_leaves{leaf}, ...
        zeros(original.npatches,1),'EdgeColor','none');
end
for leaf = 1:4
    set(surfaces(leaf),'FaceColor',colors(leaf,:),'EdgeColor','none')
end
axis(ax,'tight'); daspect(ax,[1 1 1/z_scale]); axis(ax,'off'); view(ax,35,25)
