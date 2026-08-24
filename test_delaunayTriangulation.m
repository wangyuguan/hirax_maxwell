clear
close all
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

%% Current dish geometry

dish_radius = 69.25;
thickness = 0.01*dish_radius;
surface_order = 4;

chunkie_order = 20;
chunkie_n0 = 3;
chunkie_nchs = 3;
chunkie_newton_iterations = 30;

rim_width = 0.028128271246*dish_radius;
cap_collar_width = rim_width;
cap_mesh_spacing_center = 0.2*dish_radius;
cap_mesh_spacing_side = cap_mesh_spacing_center;
cap_mesh_side_start = 0.62;
outline_refinement = 1;
wall_profile_refinement = 2;

geometry_options = struct();
geometry_options.norder = surface_order;
geometry_options.chunkie_order = chunkie_order;
geometry_options.chunkie_n0 = chunkie_n0;
geometry_options.chunkie_nchs = chunkie_nchs;
geometry_options.chunkie_newton_iterations = ...
    chunkie_newton_iterations;
geometry_options.rim_width = rim_width;
geometry_options.cap_collar_width = cap_collar_width;
geometry_options.cap_mesh_spacing_center = ...
    cap_mesh_spacing_center;
geometry_options.cap_mesh_spacing_side = ...
    cap_mesh_spacing_side;
geometry_options.cap_mesh_side_start = cap_mesh_side_start;
geometry_options.outline_refinement = outline_refinement;
geometry_options.wall_profile_refinement = ...
    wall_profile_refinement;

[~,parts] = hirax_chunkie_dish_plate_surfer( ...
    thickness,geometry_options);

%% Input points and constrained boundary edges

vertices = parts.core_vertices;
center_spacing = cap_mesh_spacing_center;
side_spacing = cap_mesh_spacing_side;

x_limits = [min(vertices(1,:)) max(vertices(1,:))];
y_limits = [min(vertices(2,:)) max(vertices(2,:))];

x_grid = x_limits(1):center_spacing:x_limits(2);
y_grid = y_limits(1):center_spacing:y_limits(2);
[grid_x,grid_y] = meshgrid(x_grid,y_grid);

inside = inpolygon(grid_x(:),grid_y(:), ...
    vertices(1,:),vertices(2,:));
interior = [grid_x(inside).';grid_y(inside).'];

distance_limit = 0.3*min(center_spacing,side_spacing);
keep = true(1,size(interior,2));
number_of_vertices = size(vertices,2);

for edge = 1:number_of_vertices
    next_edge = mod(edge,number_of_vertices)+1;
    point_a = vertices(:,edge);
    edge_vector = vertices(:,next_edge)-point_a;
    fraction = sum((interior-point_a).*edge_vector,1)/ ...
        sum(edge_vector.^2);
    fraction = min(max(fraction,0),1);
    projection = point_a+edge_vector.*fraction;
    keep = keep & ...
        vecnorm(interior-projection,2,1)>distance_limit;
end

interior = interior(:,keep);
points = [vertices interior].';
constraints = [(1:number_of_vertices).', ...
    [2:number_of_vertices 1].'];

%% Constrained Delaunay output

triangulation_object = ...
    delaunayTriangulation(points,constraints);
candidate_faces = triangulation_object.ConnectivityList;

centroids = (points(candidate_faces(:,1),:)+ ...
    points(candidate_faces(:,2),:)+ ...
    points(candidate_faces(:,3),:))/3;

use = inpolygon(centroids(:,1),centroids(:,2), ...
    vertices(1,:),vertices(2,:));
faces = candidate_faces(use,:).';
nodes = points.';

point_a = nodes(:,faces(1,:));
point_b = nodes(:,faces(2,:));
point_c = nodes(:,faces(3,:));
signed_twice_area = ...
    (point_b(1,:)-point_a(1,:)).* ...
    (point_c(2,:)-point_a(2,:))- ...
    (point_b(2,:)-point_a(2,:)).* ...
    (point_c(1,:)-point_a(1,:));

nondegenerate = abs(signed_twice_area)>1e-12;
faces = faces(:,nondegenerate);
signed_twice_area = signed_twice_area(nondegenerate);
flip = signed_twice_area<0;
faces([2 3],flip) = faces([3 2],flip);

%% Input

figure(1)
clf
hold on
closed_boundary = [vertices vertices(:,1)];
plot(closed_boundary(1,:),closed_boundary(2,:),'k-')
plot(interior(1,:),interior(2,:),'k.','MarkerSize',10)
plot(vertices(1,:),vertices(2,:),'ro', ...
    'MarkerSize',5,'LineWidth',1)
axis equal
axis tight
box on
title('Input','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%% Output

figure(2)
clf
triplot(faces.',nodes(1,:),nodes(2,:),'k')
axis equal
axis tight
box on
title('Output','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
