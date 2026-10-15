function [fig,model] = preview_dish_array(spacing,vertex_height,visible)
%PREVIEW_DISH_ARRAY Plot four zenith-pointing dishes and their four-leaf feeds.
%   [fig,model] = preview_dish_array(7,1,'on')
%
%   Dimensions are in metres. The spacing (7 m) is illustrative and the
%   vertex-ground clearance (1 m) is provisional. Dish dimensions come from
%   cad/dish_fd021_m.json. Each feed is centered at the geometric focus.
%   The leaf outline is the project's Chunkie-smoothed outline, converted
%   from mm to m and rendered as a thin extrusion for visualization.
%   Reflectors use the analytic CAD paraboloid and normal-offset back face.
%   This is a geometry illustration, not a scattering or support-frame model.
%   PNG, interactive FIG, and numeric layout data are saved in test/results.

    if nargin<1 || isempty(spacing), spacing = 7; end
    if nargin<2 || isempty(vertex_height), vertex_height = 1; end
    if nargin<3 || isempty(visible), visible = 'on'; end
    validateattributes(spacing,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(vertex_height,{'numeric'},{'scalar','real','finite','positive'});
    dish_root = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(dish_root));
    run(fullfile(project_root,'startup.m'));
    cad = jsondecode(fileread(fullfile(dish_root,'cad','dish_fd021_m.json')));
    assert(strcmp(cad.units,'m'));
    assert(spacing>2*cad.rim_outer_r,'Neighboring dishes must not overlap.');
    assert(vertex_height>cad.wall_thickness,'Back surface must clear ground.');

    model.diameter = cad.diameter;
    model.focal_ratio = cad.f_ratio;
    model.focal_length = cad.focal_length;
    model.wall_thickness = cad.wall_thickness;
    model.spacing = spacing;
    model.vertex_height = vertex_height;
    model.feed_height = vertex_height+cad.focal_length;
    model.rim_height = vertex_height+cad.diameter^2/(16*cad.focal_length);
    model.centers = spacing*[-0.5,0.5,-0.5,0.5;-0.5,-0.5,0.5,0.5];
    model.feed_centers = [model.centers;model.feed_height*ones(1,4)];
    model.leaf_thickness = 0.6925e-3;
    model.leaf_gap = 4e-3;
    model.leaf_outlines = feed_outlines(model.leaf_gap);
    model.description = ['Analytic f/D=0.21 dish and current leaf outline; ', ...
        'spacing illustrative; vertex-ground height provisional; no supports'];
    bounds = [model.leaf_outlines{:}];
    model.feed_extent = max(bounds,[],2)-min(bounds,[],2);
    assert(max(vecnorm(bounds,2,1))<cad.focal_length);

    fig = figure('Color','w','Name','HIRAX: 2 x 2 dish and four-leaf array', ...
        'NumberTitle','off','Position',[70,70,1520,850], ...
        'WindowStyle','normal','Visible',visible);
    layout = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    title(layout,'HIRAX  |  Four dishes with four-leaf feeds','FontSize',20);
    subtitle(layout,sprintf(['D = %.0f m   |   f/D = %.2f   |   f = %.2f m   |   ', ...
        'spacing %.1f m (illustrative)   |   vertex height %.1f m (provisional)'], ...
        cad.diameter,cad.f_ratio,cad.focal_length,spacing,vertex_height),'FontSize',11);

    ax = nexttile(layout,1,[2,2]);
    hold(ax,'on');
    limit = spacing/2+cad.aperture_radius+0.8;
    patch(ax,[-1,1,1,-1]*limit,[-1,-1,1,1]*limit,zeros(1,4), ...
        [0.88,0.85,0.78],'EdgeColor','none');
    grid_coordinates = -floor(limit):2:floor(limit);
    for x = grid_coordinates
      plot3(ax,[x,x],[-limit,limit],[0.002,0.002], ...
          'Color',[0.77,0.75,0.69],'LineWidth',0.5);
      plot3(ax,[-limit,limit],[x,x],[0.002,0.002], ...
          'Color',[0.77,0.75,0.69],'LineWidth',0.5);
    end
    for j = 1:4
      center = model.centers(:,j);
      draw_dish(ax,center,model);
      draw_feed(ax,model.leaf_outlines,model.leaf_thickness,model.feed_centers(:,j));
      plot3(ax,[center(1),center(1)],[center(2),center(2)], ...
          [vertex_height,model.feed_height],'--','Color',[0.77,0.32,0.14], ...
          'LineWidth',0.8);
      text(ax,center(1),center(2)-cad.aperture_radius-0.22,model.rim_height, ...
          sprintf('Dish %d',j),'HorizontalAlignment','center', ...
          'FontWeight','bold','Color',[0.15,0.25,0.31],'FontSize',11);
    end
    axis(ax,'equal'); axis(ax,'vis3d');
    xlim(ax,[-limit,limit]); ylim(ax,[-limit,limit]);
    zlim(ax,[0,model.rim_height+0.7]);
    view(ax,-36,43);
    camzoom(ax,0.75);
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'Height above ground (m)');
    title(ax,'Array overview: all components at their actual relative scale','FontSize',12);
    set(ax.Title,'Units','normalized','Position',[0.5,1.01,0]);
    set(ax,'FontSize',11,'Projection','orthographic','Box','off');
    camlight(ax,'headlight'); lighting(ax,'gouraud');
    light(ax,'Position',[-8,-6,15],'Style','infinite','Color',[0.5,0.5,0.5]);

    zoom_ax = nexttile(layout,3);
    hold(zoom_ax,'on');
    colors = [0.86,0.29,0.10;0.96,0.48,0.18;0.94,0.40,0.12;0.80,0.24,0.08];
    for j = 1:4
      xy = 100*model.leaf_outlines{j};
      patch(zoom_ax,xy(1,:),xy(2,:),colors(j,:), ...
          'EdgeColor',[0.40,0.15,0.06],'LineWidth',0.7);
    end
    plot(zoom_ax,0,0,'+','Color',[0.1,0.1,0.1],'MarkerSize',8,'LineWidth',1);
    axis(zoom_ax,'equal'); grid(zoom_ax,'on');
    extent = max(abs(bounds(:)))*100+1.5;
    xlim(zoom_ax,[-extent,extent]); ylim(zoom_ax,[-extent,extent]);
    xlabel(zoom_ax,'Local x (cm)'); ylabel(zoom_ax,'Local y (cm)');
    title(zoom_ax,{'Feed close-up: same four leaves on every dish', ...
        sprintf('%.1f x %.1f cm footprint; center at focus',100*model.feed_extent)}, ...
        'FontSize',11);
    set(zoom_ax,'FontSize',10,'Box','on');

    side_ax = nexttile(layout,6);
    hold(side_ax,'on');
    x = linspace(-cad.aperture_radius,cad.aperture_radius,301);
    z = vertex_height+x.^2/(4*cad.focal_length);
    plot(side_ax,[-3.6,3.6],[0,0],'-','Color',[0.57,0.49,0.33],'LineWidth',3);
    plot(side_ax,x,z,'Color',[0.24,0.45,0.61],'LineWidth',2.2);
    plot(side_ax,[-0.5,0.5]*model.feed_extent(1), ...
        [model.feed_height,model.feed_height],'-','Color',[0.87,0.31,0.10],'LineWidth',3);
    plot(side_ax,0,model.feed_height,'o','MarkerSize',5, ...
        'MarkerFaceColor',[0.87,0.31,0.10],'MarkerEdgeColor','none');
    plot(side_ax,[0,0],[0,model.feed_height],'--','Color',[0.5,0.5,0.5]);
    text(side_ax,0.20,vertex_height/2,sprintf('h = %.2f m',vertex_height),'FontSize',10);
    text(side_ax,0.20,vertex_height+cad.focal_length/2, ...
        sprintf('f = %.2f m',cad.focal_length),'FontSize',10);
    text(side_ax,0,model.feed_height+0.22, ...
        sprintf('Feed: z = %.2f m',model.feed_height), ...
        'HorizontalAlignment','center','Color',[0.65,0.22,0.08],'FontSize',10);
    axis(side_ax,'equal'); xlim(side_ax,[-3.6,3.6]);
    ylim(side_ax,[-0.18,model.rim_height+0.65]);
    xlabel(side_ax,'Local x (m)'); ylabel(side_ax,'z (m)');
    title(side_ax,'Single-dish section: feed is below the rim','FontSize',11);
    set(side_ax,'FontSize',10,'Box','on'); grid(side_ax,'on');

    rotate3d(fig,'on');
    drawnow;
    out_dir = fullfile(project_root,'test','results');
    if ~isfolder(out_dir), mkdir(out_dir); end
    stem = fullfile(out_dir,'hirax_array_2x2');
    model.png_file = [stem,'.png'];
    model.fig_file = [stem,'.fig'];
    model.mat_file = [stem,'.mat'];
    save(model.mat_file,'model');
    savefig(fig,model.fig_file);
    exportgraphics(fig,model.png_file,'Resolution',150);
    fprintf('Four dishes: D %.2f m, f/D %.2f, spacing %.2f m\n', ...
        model.diameter,model.focal_ratio,spacing);
    fprintf('Vertex/feed/rim heights: %.6f / %.6f / %.6f m\n', ...
        model.vertex_height,model.feed_height,model.rim_height);
    fprintf('Feed extent: %.6f x %.6f m; thickness %.7f m\n', ...
        model.feed_extent,model.leaf_thickness);
    fprintf('Saved %s\nHIRAX_ARRAY_PREVIEW_READY\n',model.fig_file);
end

function outlines = feed_outlines(gap)
    parts = leaf_master_geometry();
    t = linspace(-1,1,7);
    polynomials = lege.pols(t(1:end-1),parts.outline.order-1);
    polynomials = reshape(polynomials,parts.outline.order,[]);
    xy = zeros(2,6*parts.outline.number_of_panels);
    for j = 1:parts.outline.number_of_panels
      xy(:,6*(j-1)+(1:6)) = parts.outline.position_coefficients(:,:,j)*polynomials;
    end
    center = parts.design.leaf_center_raw;
    center(2) = center(2)+(1000*gap-2.83*sqrt(2))/sqrt(2);
    angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
    outlines = cell(1,4);
    for j = 1:4
      a = angles(j);
      outlines{j} = 1e-3*[cos(a),-sin(a);sin(a),cos(a)]*(xy+center);
    end
end

function draw_dish(ax,center,model)
    radius = linspace(0,model.diameter/2,65);
    theta = linspace(0,2*pi,145).';
    z = model.vertex_height+radius.^2/(4*model.focal_length);
    slope = radius/(2*model.focal_length);
    outer_r = radius+model.wall_thickness*slope./sqrt(1+slope.^2);
    outer_z = z-model.wall_thickness./sqrt(1+slope.^2);
    surf(ax,center(1)+cos(theta)*outer_r,center(2)+sin(theta)*outer_r, ...
        ones(size(theta))*outer_z,'FaceColor',[0.43,0.58,0.68], ...
        'EdgeColor','none','AmbientStrength',0.55,'SpecularStrength',0.25);
    surf(ax,center(1)+cos(theta)*radius,center(2)+sin(theta)*radius, ...
        ones(size(theta))*z,'FaceColor',[0.73,0.83,0.88], ...
        'EdgeColor','none','AmbientStrength',0.6,'SpecularStrength',0.25);
    surf(ax,center(1)+cos(theta)*[radius(end),outer_r(end)], ...
        center(2)+sin(theta)*[radius(end),outer_r(end)], ...
        ones(size(theta))*[z(end),outer_z(end)], ...
        'FaceColor',[0.55,0.69,0.78],'EdgeColor','none');
    plot3(ax,center(1)+radius(end)*cos(theta), ...
        center(2)+radius(end)*sin(theta),z(end)*ones(size(theta)), ...
        'Color',[0.25,0.40,0.50],'LineWidth',0.6);
end

function draw_feed(ax,outlines,thickness,center)
    for j = 1:4
      xy = outlines{j}+center(1:2);
      for sign = [-1,1]
        patch(ax,xy(1,:),xy(2,:), ...
            (center(3)+sign*thickness/2)*ones(1,size(xy,2)), ...
            [0.90,0.30,0.08],'EdgeColor','none','FaceLighting','none');
      end
      xy = xy(:,[1:end,1]);
      surf(ax,[xy(1,:);xy(1,:)],[xy(2,:);xy(2,:)], ...
          center(3)+[-1;1]*ones(1,size(xy,2))*thickness/2, ...
          'FaceColor',[0.65,0.20,0.06],'EdgeColor','none');
    end
end
