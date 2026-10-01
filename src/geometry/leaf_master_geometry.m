function parts = leaf_master_geometry(opts)
% Fixed Chunkie-smoothed master outline in the centered leaf frame.

if nargin<1 || isempty(opts), opts = struct(); end
if ~isfield(opts,'chunkie_order') || isempty(opts.chunkie_order)
    opts.chunkie_order = 20;
end
if ~isfield(opts,'chunkie_n0') || isempty(opts.chunkie_n0)
    opts.chunkie_n0 = 3;
end
if ~isfield(opts,'chunkie_nchs') || isempty(opts.chunkie_nchs)
    opts.chunkie_nchs = 3;
end
if ~isfield(opts,'chunkie_newton_iterations') || isempty(opts.chunkie_newton_iterations)
    opts.chunkie_newton_iterations = 30;
end

persistent cached_parameters cached_parts
parameters = [opts.chunkie_order opts.chunkie_n0 opts.chunkie_nchs ...
    opts.chunkie_newton_iterations];
if ~isempty(cached_parameters) && isequal(parameters,cached_parameters)
    parts = cached_parts;
    return
end

x0 = 0; y0 = 80.35;
x1 = -49.25; y1 = 80.35;
x2 = 49.25; y2 = 80.35;
x3 = 0; y3 = 2.83;
circle_radius = 20.0;
bottom_half_width = 13.8;
ellipse_radius_x = 69.25;
ellipse_radius_y = 54.41;

bottom_left = [x3-bottom_half_width;y3+bottom_half_width];
bottom_right = [x3+bottom_half_width;y3+bottom_half_width];
right_circle_start = [x2+circle_radius*cos(-pi/4); ...
    y2+circle_radius*sin(-pi/4)];
right_ellipse_join = [x2+circle_radius;y2];
left_ellipse_join = [x1-circle_radius;y1];
left_circle_end = [x1+circle_radius*cos(5*pi/4); ...
    y1+circle_radius*sin(5*pi/4)];

n0 = opts.chunkie_n0;
bottom = open_line(bottom_left,bottom_right,n0);
right_line = open_line(bottom_right,right_circle_start,ceil(7*n0/4));
right_circle = open_circle([x2;y2],circle_radius,-pi/4,0,n0);
upper_ellipse = open_ellipse([x0;y0],ellipse_radius_x, ...
    ellipse_radius_y,0,pi,ceil(30*n0/4));
left_circle = open_circle([x1;y1],circle_radius,pi,5*pi/4,n0);
left_line = open_line(left_circle_end,bottom_left,ceil(7*n0/4));
raw_vertices = [bottom right_line right_circle upper_ellipse ...
    left_circle left_line];

smoother_options.k = opts.chunkie_order;
smoother_options.n_newton = opts.chunkie_newton_iterations;
smoother_options.nchs = opts.chunkie_nchs;
master = sort(chnk.smoother.smooth(raw_vertices,smoother_options));

leaf_center = [0;0.5*(bottom_left(2)+y0+ellipse_radius_y)];
number_of_panels = master.nch;
position_coefficients = exps(master);
position_coefficients(:,1,:) = position_coefficients(:,1,:)-leaf_center;

outline.order = master.k;
outline.number_of_panels = number_of_panels;
outline.position_coefficients = position_coefficients;
outline.raw_vertices = raw_vertices-leaf_center;
outline.panel_starts = zeros(2,number_of_panels);
outline.panel_start_derivatives = zeros(2,number_of_panels);
for panel = 1:number_of_panels
    [outline.panel_starts(:,panel),outline.panel_start_derivatives(:,panel)] = ...
        chunkie_panel_values(outline,panel,0);
end

design.bottom_left = bottom_left-leaf_center;
design.bottom_right = bottom_right-leaf_center;
design.right_circle_start = right_circle_start-leaf_center;
design.right_ellipse_join = right_ellipse_join-leaf_center;
design.left_ellipse_join = left_ellipse_join-leaf_center;
design.left_circle_end = left_circle_end-leaf_center;
design.leaf_center_raw = leaf_center;
design.leaf_scale = ellipse_radius_x;

parts.outline = outline;
parts.design = design;
parts.outline_points = outline.panel_starts;
parts.outline_derivatives = outline.panel_start_derivatives;
parts.outline_position_coefficients = position_coefficients;
parts.chunkie_order = outline.order;
parts.chunkie_raw_vertices = outline.raw_vertices;
parts.leaf_radius = design.leaf_scale;

cached_parameters = parameters;
cached_parts = parts;
end

function points = open_line(point_a,point_b,number_of_points)
q = (0:number_of_points-1)/number_of_points;
points = point_a+(point_b-point_a).*q;
end

function points = open_circle(center,radius,theta_start,theta_end,number_of_points)
q = (0:number_of_points-1)/number_of_points;
theta = theta_start+(theta_end-theta_start)*q;
points = center+radius*[cos(theta);sin(theta)];
end

function points = open_ellipse(center,radius_x,radius_y,theta_start,theta_end,number_of_points)
q = (0:number_of_points-1)/number_of_points;
theta = theta_start+(theta_end-theta_start)*q;
points = center+[radius_x*cos(theta);radius_y*sin(theta)];
end

function [position,derivative] = chunkie_panel_values(outline,panel,q)
t = 2*q-1;
[polynomials,polynomial_derivatives] = lege.pols(t,outline.order-1);
polynomials = reshape(polynomials,outline.order,[]);
polynomial_derivatives = reshape(polynomial_derivatives,outline.order,[]);
n = (0:outline.order-1).';
left = t == -1;
right = t == 1;
if any(left)
    polynomial_derivatives(:,left) = ...
        ((-1).^(n+1).*n.*(n+1)/2).*ones(1,nnz(left));
end
if any(right)
    polynomial_derivatives(:,right) = ...
        (n.*(n+1)/2).*ones(1,nnz(right));
end
coefficients = outline.position_coefficients(:,:,panel);
position = coefficients*polynomials;
derivative = 2*coefficients*polynomial_derivatives;
end
