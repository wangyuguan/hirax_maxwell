function handle = plot_surfer_patch_boundaries( ...
    ax,S,line_color,line_width,patch_ids,display_offset,nedge,patch_faces)
%PLOT_SURFER_PATCH_BOUNDARIES Draw true FMM3DBIE patch edges.
%
% This draws the polynomial patch boundaries, not the auxiliary triangles
% used internally by surfer/plot for graphics rendering. If PATCH_FACES is
% supplied for triangular patches, each shared edge is drawn only once.

if nargin < 3 || isempty(line_color)
    line_color = [0.08 0.08 0.08];
end
if nargin < 4 || isempty(line_width)
    line_width = 0.55;
end
if nargin < 5 || isempty(patch_ids)
    patch_ids = 1:S.npatches;
end
if nargin < 6 || isempty(display_offset)
    display_offset = [0;0;0];
end
if nargin < 7 || isempty(nedge)
    nedge = 41;
end
if nargin < 8
    patch_faces = [];
end
display_offset = display_offset(:);
patch_ids = patch_ids(:).';

s = linspace(-1,1,nedge);
patch_types = S.iptype(patch_ids);
number_of_edges = sum(4*(patch_types == 11 | patch_types == 12)+ ...
    3*(patch_types == 1));

edge_schedule = zeros(number_of_edges,2);
edge_cursor = 1;
for patch = patch_ids
    if any(S.iptype(patch) == [11,12])
        local_edge_count = 4;
    elseif S.iptype(patch) == 1
        local_edge_count = 3;
    end
    rows = edge_cursor:edge_cursor+local_edge_count-1;
    edge_schedule(rows,:) = [repmat(patch,local_edge_count,1), ...
        (1:local_edge_count).'];
    edge_cursor = edge_cursor+local_edge_count;
end

if ~isempty(patch_faces)
    edge_vertex_columns = [1 3;1 2;3 2];
    edge_pairs = zeros(number_of_edges,2);
    for row = 1:number_of_edges
        patch = edge_schedule(row,1);
        edge = edge_schedule(row,2);
        edge_pairs(row,:) = sort(patch_faces( ...
            patch,edge_vertex_columns(edge,:)));
    end
    [~,unique_rows] = unique(edge_pairs,'rows','stable');
    edge_schedule = edge_schedule(sort(unique_rows),:);
    number_of_edges = size(edge_schedule,1);
end

curve_data = nan(3,number_of_edges*(nedge+1));
cursor = 1;

for scheduled_edge = 1:number_of_edges
    patch = edge_schedule(scheduled_edge,1);
    edge = edge_schedule(scheduled_edge,2);
    norder = S.norders(patch);
    switch S.iptype(patch)
        case {11,12}
            edges = cat(3,[s;-ones(1,nedge)], ...
                [ones(1,nedge);s],[s;ones(1,nedge)], ...
                [-ones(1,nedge);s]);
            uv = edges(:,:,edge);
            if S.iptype(patch) == 11
                basis = polytens.lege.pols(norder,uv);
            else
                basis = polytens.cheb.pols(norder,uv);
            end
            xyz = S.srccoefs{patch}(1:3,:)*basis;
        case 1
            edges = cat(3,[zeros(1,nedge);(s+1)/2], ...
                [(s+1)/2;zeros(1,nedge)], ...
                [(s+1)/2;(1-s)/2]);
            basis = koorn.pols(norder,edges(:,:,edge));
            xyz = S.srccoefs{patch}(1:3,:)*basis;
    end
    curve_data(:,cursor:cursor+nedge-1) = xyz;
    cursor = cursor+nedge+1;
end

curve_data = curve_data+display_offset;
handle = plot3(ax,curve_data(1,:),curve_data(2,:),curve_data(3,:), ...
    'Color',line_color,'LineWidth',line_width,'Clipping','off');
end
