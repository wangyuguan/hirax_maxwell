function [outlines,info] = leaf_tip_outlines(base_parts,opts)
% Keep the shared master arc and round each open tip with Chunkie.

if nargin < 2, opts = struct(); end
opts = leaf_geometry_options(opts);
base = base_parts.outline_position_coefficients;
order = size(base,2);
nbase = size(base,3);
scale = base_parts.leaf_radius/69.25;

if ~isfield(opts,'tip_attach_y') || isempty(opts.tip_attach_y)
    opts.tip_attach_y = -50*scale;
end
if ~isfield(opts,'tip_corner_width') || isempty(opts.tip_corner_width)
    opts.tip_corner_width = scale;
end
if ~isfield(opts,'tip_connector_length') || isempty(opts.tip_connector_length)
    opts.tip_connector_length = 2*scale;
end
if ~isfield(opts,'tip_rounding_eps') || isempty(opts.tip_rounding_eps)
    opts.tip_rounding_eps = 1e-11;
end

starts = base_parts.outline_points;
right = find(starts(1,:) > 0);
left = find(starts(1,:) < 0);
[~,i] = min(abs(starts(2,right)-opts.tip_attach_y));
right_panel = right(i);
[~,i] = min(abs(starts(2,left)-opts.tip_attach_y));
left_panel = left(i);
ncommon = mod(left_panel-right_panel,nbase);
common_ids = mod(right_panel-1+(0:ncommon-1),nbase)+1;
common = base(:,:,common_ids);

[right_attach,right_d] = endpoint(common(:,:,1),false);
[left_attach,left_d] = endpoint(common(:,:,end),true);
tleft = left_d/norm(left_d);
tright = right_d/norm(right_d);
lead = 2*opts.tip_connector_length;

if isfield(opts,'tip_polygons') && ~isempty(opts.tip_polygons)
    polygons = opts.tip_polygons;
else
    polygons = provisional_polygons(scale);
end

outlines = cell(1,4);
for leaf = 1:4
    vertices = [left_attach,left_attach+lead*tleft,polygons{leaf}, ...
        right_attach-lead*tright,right_attach];
    lengths = vecnorm(diff(vertices,1,2),2,1);
    widths = [0,min(opts.tip_corner_width, ...
        0.4*min(lengths(1:end-1),lengths(2:end))),0];
    cparams.ifclosed = false;
    cparams.widths = widths;
    cparams.eps = opts.tip_rounding_eps;
    pref.k = order;
    tip = chunkerpoly(vertices,cparams,pref);
    coefficients = cat(3,common,tip.exps());

    outline.order = order;
    outline.number_of_panels = size(coefficients,3);
    outline.position_coefficients = coefficients;
    outlines{leaf} = outline;
end

% Promote only bitwise-identical connector panels, without changing curves.
if opts.extend_common_outline
    nleft = 0; nright = 0;
    limit = min(cellfun(@(o)o.number_of_panels,outlines))-ncommon;
    while nleft+nright<limit && identical_panel(outlines,ncommon+nleft+1,false)
        nleft = nleft+1;
    end
    while nleft+nright<limit && identical_panel(outlines,nright+1,true)
        nright = nright+1;
    end
    for leaf = 1:4
        n = outlines{leaf}.number_of_panels;
        ids = [n-nright+1:n,1:ncommon+nleft,ncommon+nleft+1:n-nright];
        outlines{leaf}.position_coefficients = ...
            outlines{leaf}.position_coefficients(:,:,ids);
    end
    ncommon = ncommon+nleft+nright;
end
info.common_outline_panels = 1:ncommon;
end

function same = identical_panel(outlines,index,from_end)
same = true;
for leaf = 1:4
    j = index;
    if from_end, j = outlines{leaf}.number_of_panels-index+1; end
    c = outlines{leaf}.position_coefficients(:,:,j);
    if leaf == 1, reference = c; else, same = same && isequal(c,reference); end
end
end


function [r,d] = endpoint(coefficients,at_end)
n = (0:size(coefficients,2)-1).';
if at_end
    r = sum(coefficients,2);
    d = coefficients*(n.*(n+1));
else
    r = coefficients*((-1).^n);
    d = coefficients*((-1).^(n+1).*n.*(n+1));
end
end


function polygons = provisional_polygons(scale)
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
for leaf = 1:4
    polygons{leaf}(2,:) = polygons{leaf}(2,:)-75.695;
    polygons{leaf} = scale*polygons{leaf};
end
end
