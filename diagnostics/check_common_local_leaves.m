function report = check_common_local_leaves(orders,geometry_options,checks)
%CHECK_COMMON_LOCAL_LEAVES Headless geometry checks; no solve or plotting.
% report = check_common_local_leaves([4 6],geometry_options,checks)
% checks.boundary_tolerance is an absolute length used only for matching
% polynomial patch edges. This is a sampled geometric check, not a proof
% against all possible surface self-intersections.

if nargin<1 || isempty(orders), orders = [4 6]; end
if nargin<2 || isempty(geometry_options), geometry_options = struct(); end
if nargin<3 || isempty(checks), checks = struct(); end
if ~isfield(checks,'boundary_tolerance'), checks.boundary_tolerance = 1e-4; end
if ~isfield(checks,'check_patch_edges'), checks.check_patch_edges = true; end
if ~isfield(checks,'require_closed_patch_edges')
    checks.require_closed_patch_edges = true;
end
if ~isfield(checks,'flux_tolerance'), checks.flux_tolerance = 1e-4; end
validateattributes(orders,{'numeric'},{'vector','integer','positive'});
validateattributes(checks.boundary_tolerance,{'numeric'}, ...
    {'scalar','real','finite','positive'});

root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(genpath(fullfile(root,'src')))
old_visibility = get(groot,'DefaultFigureVisible');
restore_visibility = onCleanup(@()set(groot,'DefaultFigureVisible',old_visibility));
set(groot,'DefaultFigureVisible','off')
figures_before = findall(groot,'Type','figure');

thickness = 0.01*69.25;
if isfield(geometry_options,'thickness')
    thickness = geometry_options.thickness;
    geometry_options = rmfield(geometry_options,'thickness');
end
results = cell(1,numel(orders));
baseline = [];
for iorder = 1:numel(orders)
    opts = geometry_options;
    opts.norder = orders(iorder);
    timer = tic;
    [S,parts] = generate_leaves_surfer(thickness,opts);
    r = inspect_surface(S,parts,checks);
    r.order = orders(iorder);
    r.elapsed = toc(timer);
    if isempty(baseline)
        baseline = struct('outlines',{parts.outlines}, ...
            'npatches',S.npatches,'rotation',parts.rotation,'shift',parts.shift, ...
            'common_core',parts.common_core,'local_core',{parts.local_core});
        r.order_geometry_change = 0;
    else
        assert(S.npatches==baseline.npatches, ...
            'Patch topology changed when only the surface order changed.');
        assert(isequal(parts.rotation,baseline.rotation) && ...
            isequal(parts.shift,baseline.shift), ...
            'Placements changed when only the surface order changed.');
        assert(isequal(parts.common_core.nodes,baseline.common_core.nodes) && ...
            isequal(parts.common_core.faces,baseline.common_core.faces), ...
            'The common planar core mesh changed with surface order.');
        for ileaf = 1:4
            assert(isequal(parts.local_core{ileaf}.nodes,baseline.local_core{ileaf}.nodes) && ...
                isequal(parts.local_core{ileaf}.faces,baseline.local_core{ileaf}.faces), ...
                'A local planar core mesh changed with surface order.');
        end
        changes = zeros(1,4);
        for ileaf = 1:4
            a = parts.outlines{ileaf}; b = baseline.outlines{ileaf};
            assert(isequal(size(a.position_coefficients), ...
                size(b.position_coefficients)), ...
                'Outline topology changed with surface order.');
            changes(ileaf) = max(abs(a.position_coefficients- ...
                b.position_coefficients),[],'all');
        end
        r.order_geometry_change = max(changes);
        assert(r.order_geometry_change<=r.common_tolerance, ...
            'The master outline changed with surface order.');
    end
    fprintf(['ORDER %d nodes=%d patches=%d common_position_error=%.3e ' ...
        'common_normal_error=%.3e common_weight_error=%.3e\n'], ...
        r.order,r.npts,r.npatches,max(r.common_position_error), ...
        max(r.common_normal_error),max(r.common_weight_error));
    fprintf('  min_weight=%.3e min_jacobian=%.3e min_normal_alignment=%.16g\n', ...
        r.minimum_weight,r.minimum_jacobian,r.minimum_normal_alignment);
    fprintf('  common_patches_per_leaf=%d local_patches=%s common_nodes_per_leaf=%d\n', ...
        r.common_patches,mat2str(r.local_patches),r.common_nodes);
    fprintf('  local_nodes=%s common_area_fraction=%.6f common_node_fraction=%.6f\n', ...
        mat2str(r.local_nodes),r.common_area_fraction,r.common_node_fraction);
    fprintf('  leaf_areas=%s\n',mat2str(r.leaf_area,10));
    fprintf('  signed_volumes=%s relative_flux=%s\n', ...
        mat2str(r.leaf_signed_volume,10),mat2str(r.leaf_relative_flux,4));
    fprintf('  max_outline_overlap_area=%.3e order_geometry_change=%.3e\n', ...
        max(r.outline_overlap_area,[],'all'),r.order_geometry_change);
    fprintf('  minimum_sampled_interleaf_boundary_separation=%.10g\n', ...
        min(r.sampled_outline_separation,[],'all'));
    if checks.check_patch_edges
        for ileaf = 1:4
            e = r.patch_edges{ileaf};
            fprintf(['  leaf=%d edges=%d unpaired=%d nonmanifold=%d ' ...
                'collapsed=%d max_edge_gap=%.3e tolerance=%.3e\n'], ...
                ileaf,e.number_of_edges,e.unpaired,e.nonmanifold,e.collapsed, ...
                e.maximum_gap,checks.boundary_tolerance);
        end
    end
    fprintf('  elapsed=%.2fs\n',r.elapsed);
    results{iorder} = r;
    clear S parts
end
assert(isempty(setdiff(findall(groot,'Type','figure'),figures_before)), ...
    'Geometry construction unexpectedly created a figure.');
report = struct('geometry_revision',leaf_geometry_revision(), ...
    'orders',orders,'thickness',thickness,'checks',checks, ...
    'results',{results});
fprintf('COMMON_LOCAL_LEAVES_GEOMETRY_CHECKS_DONE\n');
end

function r = inspect_surface(S,parts,checks)
required = {'common_patch_ids','local_patch_ids','local_common_surfer', ...
    'local_unique_surfers','leaves','rotation','shift','outlines','outline_info'};
for k = 1:numel(required)
    assert(isfield(parts,required{k}),'Missing parts.%s.',required{k});
end
assert(numel(parts.common_patch_ids)==4 && numel(parts.local_patch_ids)==4);
assert(numel(parts.local_unique_surfers)==4 && numel(parts.leaves)==4);
assert(numel(parts.outlines)==4);
all_patches = [];
for k = 1:4
    all_patches = [all_patches,parts.common_patch_ids{k}(:).', ...
        parts.local_patch_ids{k}(:).']; %#ok<AGROW>
end
assert(isequal(sort(all_patches),1:S.npatches), ...
    'Common/local patch IDs must cover every patch exactly once.');
assert(numel(S.ixyzs)==S.npatches+1 && S.ixyzs(1)==1 && ...
    S.ixyzs(end)==S.npts+1 && all(diff(S.ixyzs)>0));
assert(all(isfinite([S.r(:);S.du(:);S.dv(:);S.n(:);S.wts(:)])), ...
    'Geometry contains nonfinite coordinates, derivatives, normals or weights.');
assert(all(S.wts>0),'Some quadrature weights are nonpositive.');
jacobian_vector = cross(S.du,S.dv,1);
jacobian = vecnorm(jacobian_vector,2,1);
assert(all(jacobian>0),'Some surface nodes have degenerate Jacobians.');
normal_norm = vecnorm(S.n,2,1);
assert(max(abs(normal_norm-1))<1e-10,'Normals are not unit length.');
alignment = sum(jacobian_vector.*S.n,1)./jacobian;
assert(min(alignment)>1-1e-10,'Normals disagree with patch orientation.');

scale = max(1,max(vecnorm(S.r,2,1)));
tolerance = 2048*eps(scale);
reference = parts.local_common_surfer;
r = struct('npts',S.npts,'npatches',S.npatches, ...
    'common_tolerance',tolerance,'minimum_weight',min(S.wts), ...
    'minimum_jacobian',min(jacobian), ...
    'minimum_normal_alignment',min(alignment));
r.common_patches = reference.npatches;
r.common_nodes = reference.npts;
r.local_patches = cellfun(@(s)s.npatches,parts.local_unique_surfers);
r.local_nodes = cellfun(@(s)s.npts,parts.local_unique_surfers);
r.common_area_fraction = 4*sum(reference.wts)/sum(S.wts);
r.common_node_fraction = 4*reference.npts/S.npts;
r.common_position_error = zeros(1,4);
r.common_normal_error = zeros(1,4);
r.common_weight_error = zeros(1,4);
r.common_derivative_error = zeros(1,4);
r.leaf_area = zeros(1,4);
r.leaf_signed_volume = zeros(1,4);
r.leaf_relative_flux = zeros(1,4);
r.patch_edges = cell(1,4);
r.outline_proper_crossings = zeros(1,4);
r.core_meshes = cell(1,5);
if all(isfield(parts,{'common_core','local_core','interface_nodes'}))
    seam = parts.interface_nodes;
    direction = seam(:,end)-seam(:,1);
    delta = seam-seam(:,1);
    assert(max(abs(direction(1)*delta(2,:)-direction(2)*delta(1,:))) ...
        <= tolerance*norm(direction),'Common/local interface must be straight.');
    for leaf = 1:4
        xy = parts.local_core{leaf}.polygon;
        region = polyshape(xy(1,:),xy(2,:),'Simplify',false);
        assert(region.NumRegions==1 && region.NumHoles==0, ...
            'Each local cap must be one complete region without holes.');
    end
    r.core_meshes{1} = inspect_core_mesh(parts.common_core, ...
        parts.interface_nodes,tolerance);
    for k = 1:4
        r.core_meshes{k+1} = inspect_core_mesh(parts.local_core{k}, ...
            parts.interface_nodes,tolerance);
        assert(all(r.core_meshes{1}.interface_edge_direction== ...
            -r.core_meshes{k+1}.interface_edge_direction), ...
            'Common and local cap meshes do not traverse their shared seam oppositely.');
    end
end
node_count = zeros(1,S.npts,'uint8');
world_polygons = cell(1,4);
world_outlines = cell(1,4);
local_polygons = cell(1,4);
for k = 1:4
    R = parts.rotation(:,:,k); shift = parts.shift(:,k);
    assert(norm(R.'*R-eye(3),'fro')<1e-12 && abs(det(R)-1)<1e-12);
    cp = parts.common_patch_ids{k}(:).';
    lp = parts.local_patch_ids{k}(:).';
    local = parts.local_unique_surfers{k};
    leaf = parts.leaves{k};
    assert(numel(cp)==reference.npatches && numel(lp)==local.npatches);
    assert(leaf.npatches==numel(cp)+numel(lp));
    ci = nodes_for_patches(S,cp); li = nodes_for_patches(S,lp);
    assert(numel(ci)==reference.npts && numel(li)==local.npts);
    assert(leaf.npts==numel(ci)+numel(li));
    node_count([ci li]) = node_count([ci li])+1;
    assert(isequal(S.iptype(cp),reference.iptype) && ...
        isequal(S.norders(cp),reference.norders));
    r.common_position_error(k) = max(vecnorm( ...
        R.'*(S.r(:,ci)-shift)-reference.r,2,1));
    r.common_normal_error(k) = max(vecnorm(R.'*S.n(:,ci)-reference.n,2,1));
    r.common_weight_error(k) = max(abs(S.wts(ci)-reference.wts));
    r.common_derivative_error(k) = max([ ...
        vecnorm(R.'*S.du(:,ci)-reference.du,2,1), ...
        vecnorm(R.'*S.dv(:,ci)-reference.dv,2,1)]);
    assert(r.common_position_error(k)<=tolerance, ...
        'Common patch positions differ after undoing the leaf rotation.');
    assert(r.common_normal_error(k)<=2048*eps, ...
        'Common patch normals differ after undoing the leaf rotation.');
    assert(r.common_weight_error(k)<=2048*eps(max(1,max(reference.wts))), ...
        'Common patch weights differ.');
    derivative_scale = max([1,vecnorm(reference.du,2,1),vecnorm(reference.dv,2,1)]);
    assert(r.common_derivative_error(k)<=2048*eps(derivative_scale), ...
        'Common patch tangent derivatives differ.');
    assert(max(vecnorm(R.'*(S.r(:,li)-shift)-local.r,2,1))<=tolerance, ...
        'Local patch IDs do not map to the saved local surface.');
    w = leaf.wts(:).';
    r.leaf_area(k) = sum(w);
    origin = mean(leaf.r,2);
    r.leaf_signed_volume(k) = sum(sum((leaf.r-origin).*leaf.n,1).*w)/3;
    r.leaf_relative_flux(k) = norm(sum(leaf.n.*w,2))/r.leaf_area(k);
    assert(r.leaf_signed_volume(k)>0,'A closed leaf has inward net orientation.');
    assert(r.leaf_relative_flux(k)<checks.flux_tolerance, ...
        'Surface-normal flux fails the closed-surface check.');

    xy = sample_outline(parts.outlines{k},9);
    r.outline_proper_crossings(k) = proper_crossings(xy,tolerance);
    assert(r.outline_proper_crossings(k)==0, ...
        'A sampled leaf outline has a proper self-intersection.');
    local_polygons{k} = polyshape(xy(1,:),xy(2,:));
    assert(local_polygons{k}.NumRegions==1 && local_polygons{k}.NumHoles==0, ...
        'A leaf footprint is disconnected or contains a hole.');
    world = R*[xy;zeros(1,size(xy,2))]+shift;
    world_outlines{k} = world(1:2,:);
    world_polygons{k} = polyshape(world(1,:),world(2,:));
    if checks.check_patch_edges
        r.patch_edges{k} = inspect_patch_edges(leaf,checks.boundary_tolerance);
        e = r.patch_edges{k};
        if checks.require_closed_patch_edges
            assert(e.unpaired==0 && e.nonmanifold==0 && e.collapsed==0, ...
                ['Leaf %d patch edges: unpaired=%d, nonmanifold=%d, ' ...
                'collapsed=%d (tolerance %.3e).'],k,e.unpaired, ...
                e.nonmanifold,e.collapsed,checks.boundary_tolerance);
            assert(e.maximum_gap<=checks.boundary_tolerance, ...
                'Leaf %d has a sampled gap along shared patch edges.',k);
        end
    end
end
assert(all(node_count==1),'Common/local node partition has omissions or overlap.');
r.outline_overlap_area = zeros(4);
r.local_outline_difference_area = zeros(4);
r.sampled_outline_separation = inf(4);
for i = 1:4
    for j = i+1:4
        overlap = area(intersect(world_polygons{i},world_polygons{j}));
        difference = area(xor(local_polygons{i},local_polygons{j}));
        r.outline_overlap_area(i,j) = overlap;
        r.local_outline_difference_area(i,j) = difference;
        separation = polyline_distance(world_outlines{i},world_outlines{j});
        r.sampled_outline_separation(i,j) = separation;
        r.sampled_outline_separation(j,i) = separation;
        assert(overlap<1e-10*scale^2, ...
            'The sampled footprints of leaves %d and %d overlap.',i,j);
        assert(difference>1e-8*scale^2, ...
            'Leaves %d and %d have indistinguishable local tip outlines.',i,j);
        assert(separation>1e-8*scale, ...
            'The sampled boundaries of leaves %d and %d touch.',i,j);
    end
end
end

function distance = polyline_distance(a,b)
% Minimum separation of sampled closed polygonal outlines. Disjointness is
% checked separately. This is not an exact closest-point curve calculation.
ae = a(:,[2:end,1]); be = b(:,[2:end,1]);
blow = min(b,be); bhigh = max(b,be); bv = be-b;
blength2 = sum(bv.^2,1);
best2 = min(sum((b-a(:,1)).^2,1));
for i = 1:size(a,2)
    low = min(a(:,i),ae(:,i)); high = max(a(:,i),ae(:,i));
    lower_distance = max(0,max(low-bhigh,blow-high));
    candidates = find(sum(lower_distance.^2,1)<=best2);
    if isempty(candidates), continue; end
    av = ae(:,i)-a(:,i); alength2 = sum(av.^2);
    start = b(:,candidates); finish = be(:,candidates); v = bv(:,candidates);
    t = max(0,min(1,sum((a(:,i)-start).*v,1)./blength2(candidates)));
    d1 = sum((a(:,i)-(start+v.*t)).^2,1);
    t = max(0,min(1,sum((ae(:,i)-start).*v,1)./blength2(candidates)));
    d2 = sum((ae(:,i)-(start+v.*t)).^2,1);
    t = max(0,min(1,sum((start-a(:,i)).*av,1)/alength2));
    d3 = sum((start-(a(:,i)+av.*t)).^2,1);
    t = max(0,min(1,sum((finish-a(:,i)).*av,1)/alength2));
    d4 = sum((finish-(a(:,i)+av.*t)).^2,1);
    best2 = min([best2,d1,d2,d3,d4]);
    if best2==0, break; end
end
distance = sqrt(best2);
end

function result = inspect_core_mesh(core,interface,tolerance)
nodes = core.nodes; faces = core.faces;
assert(size(nodes,1)==2 && size(faces,1)==3);
assert(all(faces(:)>=1 & faces(:)<=size(nodes,2)));
a = nodes(:,faces(2,:))-nodes(:,faces(1,:));
b = nodes(:,faces(3,:))-nodes(:,faces(1,:));
twice_area = a(1,:).*b(2,:)-a(2,:).*b(1,:);
assert(all(twice_area>0),'Core cap contains inverted or degenerate triangles.');
triangle_area = sum(twice_area)/2;
polygon_area = polyarea(core.polygon(1,:),core.polygon(2,:));
assert(abs(triangle_area-polygon_area)<1e-10*max(1,polygon_area), ...
    'Core triangles do not cover their polygon to roundoff accuracy.');
centers = (nodes(:,faces(1,:))+nodes(:,faces(2,:))+nodes(:,faces(3,:)))/3;
assert(all(inpolygon(centers(1,:),centers(2,:), ...
    core.polygon(1,:),core.polygon(2,:))), ...
    'A core triangle centroid lies outside its cap polygon.');
directed = [faces([1 2],:),faces([2 3],:),faces([3 1],:)].';
[keys,~,groups] = unique(sort(directed,2),'rows');
counts = accumarray(groups,1);
assert(all(counts<=2),'Nonmanifold edges occur in the planar core mesh.');
interface_ids = zeros(1,size(interface,2));
maximum_vertex_error = 0;
for k = 1:size(interface,2)
    [distance,id] = min(vecnorm(nodes-interface(:,k),2,1));
    assert(distance<=tolerance,'A cap mesh is missing a shared interface node.');
    interface_ids(k) = id;
    maximum_vertex_error = max(maximum_vertex_error,distance);
end
assert(numel(unique(interface_ids))==numel(interface_ids));
directions = zeros(1,numel(interface_ids)-1);
for k = 1:numel(directions)
    pair = interface_ids(k:k+1);
    [found,group] = ismember(sort(pair),keys,'rows');
    assert(found && counts(group)==1, ...
        'A shared interface segment is not a boundary edge of both core meshes.');
    edge = directed(find(groups==group,1),:);
    if isequal(edge,pair), directions(k) = 1; else, directions(k) = -1; end
end
result = struct('number_of_triangles',size(faces,2), ...
    'triangle_area',triangle_area,'polygon_area',polygon_area, ...
    'minimum_triangle_area',min(twice_area)/2, ...
    'interface_edge_direction',directions, ...
    'maximum_interface_vertex_error',maximum_vertex_error);
end

function ids = nodes_for_patches(S,patches)
ids = zeros(1,sum(S.ixyzs(patches+1)-S.ixyzs(patches)));
offset = 0;
for patch = patches
    current = S.ixyzs(patch):S.ixyzs(patch+1)-1;
    ids(offset+(1:numel(current))) = current;
    offset = offset+numel(current);
end
end

function xy = sample_outline(outline,nq)
q = (0:nq-1)/nq;
B = reshape(lege.pols(2*q-1,outline.order-1),outline.order,[]);
xy = reshape(pagemtimes(outline.position_coefficients,B),2,[]);
scale = max(1,max(vecnorm(xy,2,1)));
% Ignore only consecutive positions indistinguishable at floating precision.
keep = [true,vecnorm(diff(xy,1,2),2,1)>2048*eps(scale)];
xy = xy(:,keep);
if norm(xy(:,end)-xy(:,1))<=2048*eps(scale), xy(:,end) = []; end
end

function number = proper_crossings(xy,tolerance)
% Excludes adjacent edges and crossings no larger than roundoff. Polyshape
% may still simplify nearly collinear vertices without a proper crossing.
n = size(xy,2); ends = xy(:,[2:end,1]);
low = min(xy,ends); high = max(xy,ends);
number = 0;
for i = 1:n-2
    candidates = i+2:n;
    if i==1, candidates(end) = []; end
    overlap = all(high(:,candidates)>=low(:,i)-tolerance & ...
        low(:,candidates)<=high(:,i)+tolerance,1);
    candidates = candidates(overlap);
    if isempty(candidates), continue; end
    a = xy(:,i); d = ends(:,i)-a;
    b = xy(:,candidates); e = ends(:,candidates)-b;
    p = b-a; q = ends(:,candidates)-a;
    cross1 = (d(1)*p(2,:)-d(2)*p(1,:))/norm(d);
    cross2 = (d(1)*q(2,:)-d(2)*q(1,:))/norm(d);
    p = a-b; q = ends(:,i)-b;
    cross3 = (e(1,:).*p(2,:)-e(2,:).*p(1,:))./vecnorm(e,2,1);
    cross4 = (e(1,:).*q(2,:)-e(2,:).*q(1,:))./vecnorm(e,2,1);
    first = (cross1>tolerance & cross2<-tolerance) | ...
        (cross1<-tolerance & cross2>tolerance);
    second = (cross3>tolerance & cross4<-tolerance) | ...
        (cross3<-tolerance & cross4>tolerance);
    number = number+nnz(first & second);
end
end

function result = inspect_patch_edges(S,tolerance)
% Match endpoints geometrically, then compare both edge curves. Endpoint
% matching has finite tolerance because surfer stores polynomial geometry.
nq = 9;
number_of_edges = sum(3*(S.iptype==1)+4*(S.iptype~=1));
curves = zeros(3,nq,number_of_edges);
edge = 0;
[kinds,~,groups] = unique([S.norders(:),S.iptype(:)],'rows');
for g = 1:size(kinds,1)
    p = kinds(g,1); kind = kinds(g,2);
    if kind==1
        vertices = [0 1 0;0 0 1];
    else
        vertices = [-1 1 1 -1;-1 -1 1 1];
    end
    nv = size(vertices,2);
    q = linspace(0,1,nq);
    B = cell(1,nv);
    for j = 1:nv
        uv = vertices(:,j)*(1-q)+vertices(:,mod(j,nv)+1)*q;
        if kind==1
            B{j} = koorn.pols(p,uv);
        elseif kind==11
            B{j} = polytens.lege.pols(p,uv);
        elseif kind==12
            B{j} = polytens.cheb.pols(p,uv);
        else
            error('Unsupported patch type %d.',kind);
        end
    end
    for patch = find(groups==g).'
        c = S.srccoefs{patch}(1:3,:);
        for j = 1:nv
            edge = edge+1;
            curves(:,:,edge) = c*B{j};
        end
    end
end
starts = reshape(curves(:,1,:),3,[]).';
ends = reshape(curves(:,end,:),3,[]).';
[~,~,vertex] = uniquetol([starts;ends],tolerance, ...
    'ByRows',true,'DataScale',1);
keys = [vertex(1:number_of_edges),vertex(number_of_edges+1:end)];
collapsed = nnz(keys(:,1)==keys(:,2));
[~,~,group] = unique(sort(keys,2),'rows');
counts = accumarray(group,1);
paired = find(counts==2);
[~,permutation] = sort(group);
offsets = [0;cumsum(counts)];
maximum_gap = 0;
for pair = paired.'
    ids = permutation(offsets(pair)+(1:2));
    a = curves(:,:,ids(1)); b = curves(:,:,ids(2));
    if keys(ids(1),1)~=keys(ids(2),1), b = fliplr(b); end
    maximum_gap = max(maximum_gap,max(vecnorm(a-b,2,1)));
end
result = struct('number_of_edges',number_of_edges, ...
    'unpaired',nnz(counts==1),'nonmanifold',nnz(counts>2), ...
    'collapsed',collapsed,'maximum_gap',maximum_gap);
end
