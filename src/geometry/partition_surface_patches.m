function [part_id,adjacency] = partition_surface_patches(S,max_patches)
% Partition a conforming surface into groups connected through full edges.
% Each group contains at most max_patches patches. Patch corners are matched
% within 1e-8 of the bounding-box diagonal; vertex-only contacts do not count.
vertices = cat(2,S.end_pt_verts{:});
scale = norm(max(vertices,[],2)-min(vertices,[],2));
[~,~,vertex_id] = uniquetol(vertices.',1e-8*scale, ...
    'ByRows',true,'DataScale',1);
edges = zeros(size(vertices,2),2);
edge_patch = zeros(size(vertices,2),1);
offset = 0;
for patch = 1:S.npatches
    ncorner = size(S.end_pt_verts{patch},2);
    rows = offset+(1:ncorner);
    corners = vertex_id(rows);
    edges(rows,:) = [corners(:),corners([2:end,1])];
    edge_patch(rows) = patch;
    offset = offset+ncorner;
end
[~,~,edge_id] = unique(sort(edges,2),'rows');
incidence = sparse(edge_id,edge_patch,1,max(edge_id),S.npatches);
adjacency = spones(incidence.'*incidence);
adjacency(1:S.npatches+1:end) = 0;

part_id = zeros(S.npatches,1);
queue = zeros(1,min(max_patches,S.npatches));
part = 0;
while any(part_id==0)
    part = part+1;
    queue(1) = find(part_id==0,1);
    part_id(queue(1)) = part;
    head = 1;
    count = 1;
    while head<=count && count<max_patches
        neighbors = find(adjacency(queue(head),:));
        neighbors = neighbors(part_id(neighbors)==0);
        neighbors = neighbors(1:min(numel(neighbors),max_patches-count));
        part_id(neighbors) = part;
        queue(count+(1:numel(neighbors))) = neighbors;
        count = count+numel(neighbors);
        head = head+1;
    end
end
end
