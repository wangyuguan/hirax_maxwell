function [outlines,info] = hirax_local_tip_outlines(base_parts,opts)
% Shared outer arc plus four rounded tips, CCW in the centered leaf frame.
% opts.tip_polygons: four 2-by-M waypoint arrays, left attachment to right.
% Default waypoints approximate the screenshots; see provisional_polygons.

if nargin < 2 || isempty(opts)
    opts = struct();
end
base = base_parts.outline_position_coefficients;
order = size(base,2);
nbase = size(base,3);
assert(size(base,1) == 2 && order >= 6, ...
    'The master outline must have at least six Legendre coefficients.')
if isfield(opts,'chunkie_order') && ~isempty(opts.chunkie_order)
    assert(opts.chunkie_order == order, ...
        'Keep chunkie_order equal to the master order to preserve common panels.')
end
scale = 1;
if isfield(base_parts,'leaf_radius')
    scale = base_parts.leaf_radius/69.25;
end
attach_y = option(opts,'tip_attach_y',-50*scale);
corner_width = option(opts,'tip_corner_width',1.0*scale);
connector_length = option(opts,'tip_connector_length',2*scale);
rounding_eps = option(opts,'tip_rounding_eps',1e-11);
validateattributes(corner_width,{'numeric'},{'scalar','positive','finite'})
validateattributes(connector_length,{'numeric'},{'scalar','positive','finite'})
validateattributes(rounding_eps,{'numeric'},{'scalar','positive','finite'})

starts = zeros(2,nbase);
for j = 1:nbase
    starts(:,j) = panel_jets(base(:,:,j),0);
end
right_candidates = find(starts(1,:) > 0);
left_candidates = find(starts(1,:) < 0);
assert(~isempty(right_candidates) && ~isempty(left_candidates), ...
    'The master outline must cross both sides of the local y axis.')
[~,ir] = min(abs(starts(2,right_candidates)-attach_y));
[~,il] = min(abs(starts(2,left_candidates)-attach_y));
right_panel = right_candidates(ir);
left_panel = left_candidates(il);
ncommon = mod(left_panel-right_panel,nbase);
assert(ncommon > 0 && ncommon < nbase, ...
    'The selected attachment points must bound a nonempty common arc.')
common_ids = mod(right_panel-1+(0:ncommon-1),nbase)+1;
common = base(:,:,common_ids);
[right_attach,right_d,right_d2] = panel_jets(common(:,:,1),0);
[left_attach,left_d,left_d2] = panel_jets(common(:,:,end),1);
assert(left_attach(1) < right_attach(1), ...
    'The common arc must run from the right attachment to the left attachment.')
tleft = left_d/norm(left_d);
tright = right_d/norm(right_d);
left_curvature = normal_acceleration(left_d,left_d2);
right_curvature = normal_acceleration(right_d,right_d2);
left_inner = left_attach+connector_length*tleft;
right_inner = right_attach-connector_length*tright;

if isfield(opts,'tip_polygons') && ~isempty(opts.tip_polygons)
    polygons = opts.tip_polygons;
    assert(iscell(polygons) && numel(polygons) == 4, ...
        'tip_polygons must be a four-element cell array.')
    polygons = reshape(polygons,1,4);
    approximate_defaults = false;
else
    polygons = provisional_polygons(scale);
    approximate_defaults = true;
end

[q,~,v2c] = lege.exps(order);
q = (q(:).'+1)/2;
left_connector = quintic_connector( ...
    left_attach,left_inner,connector_length*tleft, ...
    connector_length*tleft,connector_length^2*left_curvature,[0;0],q);
right_connector = quintic_connector( ...
    right_inner,right_attach,connector_length*tright, ...
    connector_length*tright,[0;0],connector_length^2*right_curvature,q);
left_coeff = (v2c*left_connector.').';
right_coeff = (v2c*right_connector.').';

outlines = cell(1,4);
info = struct();
info.common_outline_panels = 1:ncommon;
info.common_master_panel_ids = common_ids;
info.local_outline_panels = cell(1,4);
info.attachment_points = [left_attach,right_attach];
info.attachment_y_requested = attach_y;
info.attachment_tangents = [tleft,tright];
info.tip_polygons = polygons;
info.approximate_image_informed_defaults = approximate_defaults;
info.geometry_source = 'caller-supplied local polygon waypoints';
if approximate_defaults
    info.geometry_source = ...
        'provisional image-informed tips; not dimensionally digitized Gerber';
end
info.smoother = 'Chunkie chunkerpoly local Gaussian corner rounding';
info.corner_width = corner_width;
info.connector_length = connector_length;
info.join_position_error = zeros(1,4);
info.join_unit_tangent_error = zeros(1,4);
info.splice_position_error = zeros(1,4);
info.splice_unit_tangent_error = zeros(1,4);
info.common_coefficient_error = zeros(1,4);
info.number_of_panels = zeros(1,4);
info.minimum_local_curvature_radius = zeros(1,4);

for leaf = 1:4
    waypoints = polygons{leaf};
    validateattributes(waypoints,{'numeric'},{'2d','finite','nonempty'})
    assert(size(waypoints,1) == 2, ...
        'Every local tip polygon must have two coordinate rows.')
    % Extra collinear endpoint segments make the open corner-rounder's end
    % tangent agree with each connector, regardless of the first waypoint.
    vertices = [left_inner,left_inner+connector_length*tleft, ...
        waypoints,right_inner-connector_length*tright,right_inner];
    lengths = vecnorm(diff(vertices,1,2),2,1);
    assert(all(lengths > 1e-10*scale), ...
        'Consecutive local polygon vertices must be distinct.')
    widths = [0,min(corner_width,0.4*min( ...
        lengths(1:end-1),lengths(2:end))),0];
    cparams = struct('ifclosed',false,'rounded',true, ...
        'widths',widths,'eps',rounding_eps);
    pref = struct('k',order,'dim',2);
    local_curve = chunkerpoly(vertices,cparams,pref);
    nlocal = local_curve.nch+2;
    local_coeff = zeros(2,order,nlocal);
    local_coeff(:,:,1) = left_coeff;
    for j = 1:local_curve.nch
        local_coeff(:,:,j+1) = ...
            (v2c*squeeze(local_curve.r(:,:,j)).').';
    end
    local_coeff(:,:,end) = right_coeff;
    coefficients = cat(3,common,local_coeff);
    outline = struct('order',order, ...
        'number_of_panels',ncommon+nlocal, ...
        'position_coefficients',coefficients, ...
        'raw_vertices',[starts(:,common_ids),left_attach,vertices]);
    [outline.panel_starts,outline.panel_start_derivatives, ...
        join_gap,join_tangent] = joins(coefficients);
    outlines{leaf} = outline;
    info.local_outline_panels{leaf} = ncommon+(1:nlocal);
    info.join_position_error(leaf) = max(join_gap);
    info.join_unit_tangent_error(leaf) = max(join_tangent);
    info.splice_position_error(leaf) = max(join_gap([ncommon,end]));
    info.splice_unit_tangent_error(leaf) = max(join_tangent([ncommon,end]));
    delta = coefficients(:,:,1:ncommon)-common;
    info.common_coefficient_error(leaf) = max(abs(delta(:)));
    info.number_of_panels(leaf) = outline.number_of_panels;
    info.minimum_local_curvature_radius(leaf) = minimum_radius(local_coeff);
end
geometry_scale = max(1,max(vecnorm(starts,2,1)));
assert(max(info.splice_position_error) < 1e-9*geometry_scale, ...
    'A local tip did not meet its common-arc attachment point.')
assert(max(info.splice_unit_tangent_error) < 1e-7, ...
    'A local tip did not meet its common-arc tangent direction.')
end


function value = option(opts,name,default)
value = default;
if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
end
end


function [r,d,d2] = panel_jets(coefficients,q)
% Exact endpoint formulas for Legendre polynomials and their first two jets.
n = (0:size(coefficients,2)-1).';
if q == 0
    p = (-1).^n;
    dp = (-1).^(n+1).*n.*(n+1)/2;
    d2p = (-1).^n.*(n-1).*n.*(n+1).*(n+2)/8;
else
    p = ones(size(n));
    dp = n.*(n+1)/2;
    d2p = (n-1).*n.*(n+1).*(n+2)/8;
end
r = coefficients*p;
d = 2*coefficients*dp;
d2 = 4*coefficients*d2p;
end


function acceleration = normal_acceleration(d,d2)
speed = norm(d);
assert(speed > 0,'The master outline has a zero endpoint tangent.')
tangent = d/speed;
acceleration = (d2-tangent*dot(tangent,d2))/speed^2;
end


function r = quintic_connector(p0,p1,d0,d1,a0,a1,q)
% Endpoint value, tangent, and curvature matched quintic in q in [0,1].
c0 = p0;
c1 = d0;
c2 = a0/2;
u = p1-c0-c1-c2;
v = d1-c1-2*c2;
w = a1-2*c2;
c3 = 10*u-4*v+w/2;
c4 = -15*u+7*v-w;
c5 = 6*u-3*v+w/2;
r = c0+c1*q+c2*q.^2+c3*q.^3+c4*q.^4+c5*q.^5;
end


function [starts,derivatives,gaps,tangent_errors] = joins(coefficients)
npan = size(coefficients,3);
starts = zeros(2,npan);
derivatives = starts;
ends = starts;
end_derivatives = starts;
for j = 1:npan
    [starts(:,j),derivatives(:,j)] = panel_jets(coefficients(:,:,j),0);
    [ends(:,j),end_derivatives(:,j)] = panel_jets(coefficients(:,:,j),1);
end
next = [2:npan,1];
gaps = vecnorm(ends-starts(:,next),2,1);
tangent_errors = vecnorm(end_derivatives./vecnorm(end_derivatives,2,1)- ...
    derivatives(:,next)./vecnorm(derivatives(:,next),2,1),2,1);
end


function polygons = provisional_polygons(scale)
% Provisional connected outer contours traced from the supplied top view.
% Interior nested rectangles in the screenshots are not extra boundaries.
polygons = cell(1,4);
polygons{1} = [-12.2,15.0;-12.2,4.85;-8.9,4.85; ...
    -8.9,13.9;11.1,13.9].';
polygons{2} = [-11.1,13.9;8.9,13.9;8.9,4.85; ...
    12.2,4.85;12.2,15.0].';
polygons{3} = [-12.2,15.0;-12.2,4.85;-8.8,4.85; ...
    -8.8,16.4;-4.2,16.4;-4.2,13.9;4.9,13.9;4.9,16.4; ...
    8.8,16.4;8.8,13.9;11.1,13.9].';
polygons{4} = [-11.1,13.9;-8.8,13.9;-8.8,16.3; ...
    -4.8,16.3;-4.8,13.9;4.5,13.9;4.5,16.4;8.9,16.4; ...
    8.9,13.9;11.1,13.9].';
for j = 1:4
    polygons{j}(2,:) = polygons{j}(2,:)-75.695;
    polygons{j} = scale*polygons{j};
end
end


function radius = minimum_radius(coefficients)
k = size(coefficients,2);
t = linspace(-1,1,max(101,4*k));
p = zeros(k,numel(t));
dp = p;
d2p = p;
p(1,:) = 1;
p(2,:) = t;
dp(2,:) = 1;
for n = 2:k-1
    p(n+1,:) = ((2*n-1)*t.*p(n,:)-(n-1)*p(n-1,:))/n;
    dp(n+1,:) = ((2*n-1)*(p(n,:)+t.*dp(n,:))- ...
        (n-1)*dp(n-1,:))/n;
    d2p(n+1,:) = ((2*n-1)*(2*dp(n,:)+t.*d2p(n,:))- ...
        (n-1)*d2p(n-1,:))/n;
end
max_curvature = 0;
for j = 1:size(coefficients,3)
    d = coefficients(:,:,j)*dp;
    d2 = coefficients(:,:,j)*d2p;
    curvature = abs(d(1,:).*d2(2,:)-d(2,:).*d2(1,:))./ ...
        vecnorm(d,2,1).^3;
    max_curvature = max(max_curvature,max(curvature));
end
radius = 1/max_curvature;
end
