%% Wall-profile data

clear
close all
clc

% Physical dimensions used by the HIRAX leaf.
leaf_radius = 69.25;
thickness = 0.01*leaf_radius;
rim_width = 0.028128271246*leaf_radius;
number_of_samples = 401;

% Quintic Bezier control points for the upper wall profile.
radial_control = [0 0.4 0.8 1 1 1];
height_control = [1 1 1 0.8 0.4 0];
control_points = [radial_control;height_control];

s = linspace(0,1,number_of_samples);
normalized_wall = evaluate_bezier(control_points,s);

half_thickness = thickness/2;
upper_x = rim_width*normalized_wall(1,:);
upper_z = half_thickness*normalized_wall(2,:);
lower_x = upper_x;
lower_z = -upper_z;
physical_control_points = [rim_width*radial_control; ...
    half_thickness*height_control];

%% Physical wall cross-section

wall_figure = figure('Color','w', ...
    'Name','HIRAX wall and Bezier profile', ...
    'Position',[80 120 960 440]);

% Physical top-cap, wall, bottom-cap, and Bezier control points.
cap_line_length = 0.85*rim_width;
axes_one = axes(wall_figure);
hold(axes_one,'on')

plot(axes_one,[-cap_line_length 0], ...
    [half_thickness half_thickness],'k-','LineWidth',2.0)
plot(axes_one,physical_control_points(1,:), ...
    physical_control_points(2,:),'k--','LineWidth',1.2)
plot(axes_one,upper_x,upper_z,'b-','LineWidth',2.4)
plot(axes_one,lower_x,lower_z,'b-','LineWidth',2.4)
plot(axes_one,[-cap_line_length 0], ...
    [-half_thickness -half_thickness],'k-','LineWidth',2.0)

scatter(axes_one,[0 0 rim_width], ...
    [half_thickness -half_thickness 0],48,'k','filled')
scatter(axes_one,physical_control_points(1,:), ...
    physical_control_points(2,:),70,'o', ...
    'MarkerEdgeColor','r','LineWidth',1.2)
for point = 1:6
    text(axes_one,physical_control_points(1,point)+0.025*rim_width, ...
        physical_control_points(2,point)+0.06*half_thickness, ...
        sprintf('$P_{%d}$',point-1), ...
        'Color','r','FontSize',11,'Interpreter','latex')
end

axis(axes_one,'equal')
xlim(axes_one,[-1.02*cap_line_length 1.63*rim_width])
ylim(axes_one,1.85*[-half_thickness half_thickness])
axis(axes_one,'off')

fprintf('Thickness: %.12g\n',thickness)
fprintf('Half thickness: %.12g\n',half_thickness)
fprintf('Rim width: %.12g\n',rim_width)


%% Local function

function curve = evaluate_bezier(control_points,s)
degree = size(control_points,2)-1;
curve = zeros(size(control_points,1),numel(s));
for k = 0:degree
    basis = nchoosek(degree,k)*(1-s).^(degree-k).*s.^k;
    curve = curve+control_points(:,k+1)*basis;
end
end
