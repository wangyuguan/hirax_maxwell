function [fig,model] = preview_toy_antenna_array(visible)
%PREVIEW_TOY_ANTENNA_ARRAY Show a low-cost surrogate for the dish/feed array.
%   [fig,model] = preview_toy_antenna_array('on')
%   Only analytic ellipsoids are sampled for drawing; no solver or MEX is used.
%   See toy_antenna_array_geometry(4) to build a coarse BIE-ready station.

    if nargin<1 || isempty(visible), visible = 'on'; end
    model = toy_antenna_array_geometry(0);
    fig = figure('Name','Toy array: closed bodies and nested interfaces', ...
        'NumberTitle','off','Color','w','WindowStyle','normal', ...
        'Position',[80,80,1500,840],'Visible',visible);
    caption(fig,[0.02,0.94,0.96,0.045], ...
        'Reduced array model  |  Five closed surfaces per antenna',19);
    caption(fig,[0.02,0.90,0.96,0.035], ...
        ['h / D = 1 / 6: equivalent to 1 m clearance for a 6 m dish; ', ...
        'ground plane is a reference only'],11);

    ax = axes(fig,'Position',[0.06,0.18,0.56,0.63]); hold(ax,'on');
    limit = model.spacing/2+max(model.outer_axes(1:2))+0.2;
    ground = patch(ax,[-1,1,1,-1]*limit,[-1,-1,1,1]*limit,model.ground_z*ones(1,4), ...
        [0.90,0.88,0.83],'EdgeColor','none','FaceLighting','none');
    for j = 1:4
      shift = model.station_shifts(:,j);
      draw_station(ax,model,shift,true);
      text(ax,shift(1),shift(2)-model.outer_axes(2)-0.08,0.06, ...
          sprintf('Antenna %d',j),'HorizontalAlignment','center', ...
          'FontSize',11,'Color',[0.17,0.25,0.30]);
    end
    axis(ax,'equal'); axis(ax,'vis3d');
    ztop = model.outer_center(3)+model.outer_axes(3)+0.1;
    xlim(ax,[-limit,limit]); ylim(ax,[-limit,limit]); zlim(ax,[model.ground_z,ztop]);
    view(ax,-35,28); camzoom(ax,1.0);
    xlabel(ax,'x / wavelength'); ylabel(ax,'y / wavelength'); zlabel(ax,'z / wavelength');
    caption(fig,[0.05,0.83,0.58,0.04],'2 x 2 array: independent outer interfaces',12);
    set(ax,'FontSize',11,'Projection','orthographic');
    camlight(ax,'headlight'); lighting(ax,'gouraud');
    set(ground,'FaceLighting','none');

    single_ax = axes(fig,'Position',[0.72,0.52,0.24,0.29]); hold(single_ax,'on');
    draw_station(single_ax,model,zeros(3,1),false);
    % Only meridian outlines in the side view: edge-on dashed circles can
    % produce export artifacts when their projected area is exactly zero.
    t = linspace(0,2*pi,161);
    plot3(single_ax,model.outer_axes(1)*cos(t),zeros(size(t)), ...
        model.outer_center(3)+model.outer_axes(3)*sin(t), ...
        '--','Color',[0.18,0.52,0.36],'LineWidth',0.9);
    plot3(single_ax,model.inner_radius*cos(t),zeros(size(t)), ...
        model.inner_center(3)+model.inner_radius*sin(t), ...
        '--','Color',[0.71,0.38,0.08],'LineWidth',0.9);
    plot3(single_ax,[-0.95,0.95],[0,0],model.ground_z*[1,1], ...
        '-','Color',[0.45,0.42,0.37],'LineWidth',1.5);
    bottom = model.body_center(3)-model.body_axes(3);
    plot3(single_ax,[0,0],[0,0],[model.ground_z,bottom], ...
        '-','Color',[0.15,0.15,0.15],'LineWidth',1.2);
    plot3(single_ax,[-0.04,0.04,nan,-0.04,0.04],zeros(1,5), ...
        [model.ground_z,model.ground_z,nan,bottom,bottom], ...
        '-','Color',[0.15,0.15,0.15],'LineWidth',1.2);
    text(single_ax,0.08,0,(model.ground_z+bottom)/2,'h = D/6', ...
        'FontSize',10,'HorizontalAlignment','left');
    axis(single_ax,'equal'); axis(single_ax,'vis3d');
    xlim(single_ax,[-0.95,0.95]); ylim(single_ax,[-0.9,0.9]); zlim(single_ax,[-0.04,ztop]);
    view(single_ax,0,0);
    caption(fig,[0.68,0.82,0.31,0.065],{'Side view: ground clearance', ...
        'Outer ellipsoid stays above ground'},11);
    axis(single_ax,'off');
    camlight(single_ax,'headlight'); lighting(single_ax,'gouraud');

    leaf_ax = axes(fig,'Position',[0.72,0.13,0.24,0.28]); hold(leaf_ax,'on');
    for j = 1:4
      a = model.leaf_angles(j);
      R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
      draw_ellipsoid(leaf_ax,model.leaf_axes, ...
          model.leaf_centers(:,j)-model.inner_center,R,[0.91,0.35,0.12]);
    end
    t = linspace(0,2*pi,160);
    plot3(leaf_ax,model.inner_radius*cos(t),model.inner_radius*sin(t),zeros(size(t)), ...
        '--','Color',[0.71,0.38,0.08],'LineWidth',1.2);
    axis(leaf_ax,'equal'); view(leaf_ax,0,90);
    xlim(leaf_ax,[-0.32,0.32]); ylim(leaf_ax,[-0.32,0.32]);
    caption(fig,[0.68,0.43,0.31,0.065],{'Four smooth closed ellipsoids', ...
        'Diameter 0.17; enclosing radius 0.28'},11);
    xlabel(leaf_ax,'Local x / wavelength'); ylabel(leaf_ax,'Local y / wavelength');
    set(leaf_ax,'FontSize',10); grid(leaf_ax,'on');
    camlight(leaf_ax,'headlight'); lighting(leaf_ax,'gouraud');

    % Dummy handles give a single legend without picking up sphere circles.
    h1 = plot3(ax,nan,nan,nan,'s','MarkerFaceColor',[0.34,0.60,0.76], ...
        'MarkerEdgeColor','none','MarkerSize',9);
    h2 = plot3(ax,nan,nan,nan,'s','MarkerFaceColor',[0.91,0.35,0.12], ...
        'MarkerEdgeColor','none','MarkerSize',9);
    h3 = plot3(ax,nan,nan,nan,'--','Color',[0.71,0.38,0.08],'LineWidth',1.2);
    h4 = plot3(ax,nan,nan,nan,'--','Color',[0.18,0.52,0.36],'LineWidth',1.0);
    key = legend(ax,[h1,h2,h3,h4],{'Body surrogate','Leaf surrogates', ...
        'Feed interface','Whole-antenna interface'}, ...
        'NumColumns',2,'FontSize',10);
    set(key,'Units','normalized','Position',[0.16,0.035,0.38,0.06]);
    rotate3d(fig,'on'); drawnow;
    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    out_dir = fullfile(project_root,'test','results');
    if ~isfolder(out_dir), mkdir(out_dir); end
    stem = fullfile(out_dir,'toy_antenna_array');
    model.png_file = [stem,'.png']; model.fig_file = [stem,'.fig'];
    save([stem,'.mat'],'model');
    savefig(fig,model.fig_file);
    exportgraphics(fig,model.png_file,'Resolution',150);
    fprintf('Saved %s\nTOY_ANTENNA_ARRAY_PREVIEW_READY\n',model.fig_file);
end

function caption(fig,position,label,size)
    annotation(fig,'textbox',position,'String',label,'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',size,'FontWeight','bold','Interpreter','none');
end

function draw_station(ax,model,shift,show_proxies)
    draw_ellipsoid(ax,model.body_axes,model.body_center+shift,eye(3),[0.34,0.60,0.76]);
    for j = 1:4
      a = model.leaf_angles(j);
      R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
      draw_ellipsoid(ax,model.leaf_axes,model.leaf_centers(:,j)+shift,R,[0.91,0.35,0.12]);
    end
    if show_proxies
      draw_interface(ax,model.inner_center+shift,model.inner_radius*ones(3,1),[0.71,0.38,0.08]);
      draw_interface(ax,model.outer_center+shift,model.outer_axes,[0.18,0.52,0.36]);
    end
end

function draw_ellipsoid(ax,abc,center,R,color)
    [x,y,z] = sphere(40);
    xyz = R*(abc.*[x(:).';y(:).';z(:).'])+center;
    surf(ax,reshape(xyz(1,:),size(x)),reshape(xyz(2,:),size(y)), ...
        reshape(xyz(3,:),size(z)),'FaceColor',color,'EdgeColor','none', ...
        'AmbientStrength',0.55,'SpecularStrength',0.22);
end

function draw_interface(ax,center,abc,color)
    t = linspace(0,2*pi,161);
    for plane = 1:3
      xyz = zeros(3,numel(t));
      coords = setdiff(1:3,plane);
      xyz(coords,:) = abc(coords).*[cos(t);sin(t)];
      xyz = xyz+center;
      plot3(ax,xyz(1,:),xyz(2,:),xyz(3,:),'--','Color',color,'LineWidth',0.9);
    end
end
