function [J,blocks,Pmerge,info] = multilevel_nrccie_skeleton(S,samplefun,tol,max_patches,direction)
%MULTILEVEL_NRCCIE_SKELETON Complete-node IDs on a spatial binary merge tree.
%
% Syntax:
%   [J,blocks,Pmerge,info] = multilevel_nrccie_skeleton( ...
%       S,samplefun,tol,max_patches,direction)
%
% Input arguments:
%   S           - source surfer, including patch connectivity and weights.
%   samplefun   - @(T) returning actual samples on a subset of surface nodes:
%                 3*T.npts rows for incoming, 3*T.npts columns for outgoing.
%                 The other dimension denotes fixed samples at all levels.
%   tol         - relative weighted Frobenius sample-error tolerance.
%   max_patches - maximum patches per connected finest-level region.
%   direction   - 'incoming' (row ID) or 'outgoing' (column ID).
%
% Output arguments:
%   J      - final skeleton DOFs, three consecutive DOFs per selected node.
%   blocks - local nodes, interpolation P, and candidate_dofs per region.
%   Pmerge - composed map: final values to local candidates for incoming;
%            local candidate densities to final skeletons for outgoing.
%   info   - sample-error bound, hierarchy, storage estimates and timings.
%
% Spatial binary splits group connected regions; each deterministic QR merge
% samples only its children's skeleton nodes. Per-level error budgets account
% for child interpolation amplification using weighted Gram-matrix bounds.
% These bounds concern supplied samples, not unseen incident or exterior data.
% Pmerge is generally dense; high-rank regions can retain all candidates.

validateattributes(tol,{'numeric'},{'real','scalar','finite','>',0,'<',1})
validateattributes(max_patches,{'numeric'},{'real','scalar','integer','positive'})
assert(any(strcmp(direction,{'incoming','outgoing'})))
assert(all(S.wts>0),'Skeletonization requires positive source weights.')
incoming = strcmp(direction,'incoming');
timer = tic;
[part_id,adjacency] = partition_surface_patches(S,max_patches);
nparts = max(part_id);
assert(nparts>0,'The surface must contain at least one patch.')
node_part = part_id(S.patch_id);
nlevels = ceil(log2(nparts))+1;
level_tol = tol/nlevels;
blocks = repmat(struct('nodes',[],'P',[],'skeleton_nodes',[], ...
    'candidate_dofs',[]),nparts,1);
part_node_counts = zeros(nparts,1);
local_skeleton_counts = zeros(nparts,1);
local_errors = zeros(nparts,1);
centers = zeros(3,nparts);
candidate_nodes = [];
leaves = cell(1,nparts);
tree = struct('children',{},'level',{},'part_count',{},'input_nodes',{}, ...
    'skeleton_nodes',{},'requested_tol',{},'local_error',{}, ...
    'error_norm_bound',{},'interpolation_norm_bound',{},'seconds',{});
largest_sample_nodes = 0;
sample_count = [];
sample_seconds = 0;
id_seconds = 0;
interpolation_entries = 0;

for part = 1:nparts
    leaf_timer = tic;
    nodes = find(node_part==part).';
    centers(:,part) = mean(S.r(:,nodes),2);
    [source,A,sample_norm] = samples(nodes);
    [selected,P,local_error] = local_id(source,A,level_tol);
    selected_nodes = nodes((selected(1:3:end)+2)/3);
    offsets = 3*numel(candidate_nodes)+(1:numel(selected));
    blocks(part).nodes = nodes;
    blocks(part).P = P;
    blocks(part).skeleton_nodes = selected_nodes;
    blocks(part).candidate_dofs = offsets;
    candidate_nodes = [candidate_nodes,selected_nodes];
    part_node_counts(part) = numel(nodes);
    local_skeleton_counts(part) = numel(selected_nodes);
    local_errors(part) = local_error;
    interpolation_entries = interpolation_entries+numel(P);
    L = weighted_map(source,P,selected);
    if incoming, gram = L'*L; else, gram = L*L'; end
    leaf = struct('nodes',selected_nodes,'map',speye(numel(selected)), ...
        'origin_dofs',offsets,'gram',gram,'sample_norm',sample_norm, ...
        'error_norm',local_error*sample_norm,'level',0,'id',part);
    leaves{part} = leaf;
    tree(part) = node_record([],0,1,numel(nodes),numel(selected_nodes), ...
        level_tol,local_error,leaf.error_norm,gram,toc(leaf_timer));
    clear A P L gram leaf
end

root = merge_region(1:nparts);
[~,order] = sort(root.origin_dofs);
if incoming
    Pmerge = root.map(order,:);
else
    Pmerge = root.map(:,order);
end
J = reshape(3*root.nodes+(-2:0).',1,[]);
relative_error_bound = root.error_norm/max(root.sample_norm,realmin);
assert(relative_error_bound<=1.1*tol, ...
    'The composed multilevel sample error bound exceeds the tolerance.')
info = struct('part_id',part_id,'adjacency',adjacency, ...
    'part_patch_counts',accumarray(part_id,1), ...
    'part_node_counts',part_node_counts, ...
    'local_skeleton_counts',local_skeleton_counts, ...
    'candidate_nodes',candidate_nodes,'local_errors',local_errors, ...
    'relative_error_bound',relative_error_bound,'tree',tree, ...
    'levels',root.level+1,'randomized',false, ...
    'largest_sample_nodes',largest_sample_nodes, ...
    'full_sample_mib',16*(3*S.npts)*sample_count/2^20, ...
    'largest_sample_mib',16*(3*largest_sample_nodes)*sample_count/2^20, ...
    'interpolation_mib',16*(interpolation_entries+numel(Pmerge))/2^20, ...
    'sample_seconds',sample_seconds,'id_seconds',id_seconds, ...
    'seconds',toc(timer));

    function node = merge_region(parts)
        if numel(parts)==1
            node = leaves{parts};
            return
        end
        span = max(centers(:,parts),[],2)-min(centers(:,parts),[],2);
        [~,axis] = max(span);
        % Part index is an explicit deterministic tie-breaker.
        [~,order_parts] = sortrows([centers(axis,parts).',parts.'],[1,2]);
        parts = parts(order_parts);
        split = floor(numel(parts)/2);
        left = merge_region(parts(1:split));
        right = merge_region(parts(split+1:end));
        merge_timer = tic;
        nodes = [left.nodes,right.nodes];
        region_norm = hypot(left.sample_norm,right.sample_norm);
        child_error = hypot(left.error_norm,right.error_norm);
        amplification = max(gram_bound(left.gram),gram_bound(right.gram));
        if isempty(nodes)
            selected = [];
            P = zeros(0,0);
            L = P;
            sample_norm = 0;
            local_error = 0;
            merge_tol = level_tol;
        else
            [source,A,sample_norm] = samples(nodes);
            merge_tol = min(level_tol,level_tol*region_norm/ ...
                max(amplification*sample_norm,realmin));
            [selected,P,local_error] = local_id(source,A,merge_tol);
            L = weighted_map(source,P,selected);
            clear A
        end
        left_dofs = 1:3*numel(left.nodes);
        right_dofs = 3*numel(left.nodes)+(1:3*numel(right.nodes));
        if incoming
            map = [left.map*P(left_dofs,:);right.map*P(right_dofs,:)];
            gram = L(left_dofs,:)'*left.gram*L(left_dofs,:) ...
                +L(right_dofs,:)'*right.gram*L(right_dofs,:);
        else
            map = [P(:,left_dofs)*left.map,P(:,right_dofs)*right.map];
            gram = L(:,left_dofs)*left.gram*L(:,left_dofs)' ...
                +L(:,right_dofs)*right.gram*L(:,right_dofs)';
        end
        selected_nodes = nodes((selected(1:3:end)+2)/3);
        error_norm = child_error+amplification*local_error*sample_norm;
        level = max(left.level,right.level)+1;
        id = numel(tree)+1;
        tree(id) = node_record([left.id,right.id],level,numel(parts), ...
            numel(nodes),numel(selected_nodes),merge_tol,local_error, ...
            error_norm,gram,toc(merge_timer));
        node = struct('nodes',selected_nodes,'map',map, ...
            'origin_dofs',[left.origin_dofs,right.origin_dofs], ...
            'gram',gram,'sample_norm',region_norm,'error_norm',error_norm, ...
            'level',level,'id',id);
    end

    function [source,A,sample_norm] = samples(nodes)
        source = node_subset(S,nodes);
        sample_timer = tic;
        A = samplefun(source);
        sample_seconds = sample_seconds+toc(sample_timer);
        w = sqrt(repelem(source.wts(:),3));
        if incoming
            assert(size(A,1)==3*numel(nodes),'Incoming sample row count mismatch.')
            count = size(A,2);
            sample_norm = norm(w.*A,'fro');
        else
            assert(size(A,2)==3*numel(nodes),'Outgoing sample column count mismatch.')
            count = size(A,1);
            sample_norm = norm(A./w.','fro');
        end
        if isempty(sample_count), sample_count = count; end
        assert(count==sample_count,'All levels must use the same sample columns/rows.')
        largest_sample_nodes = max(largest_sample_nodes,numel(nodes));
    end

    function [selected,P,error] = local_id(source,A,local_tol)
        id_timer = tic;
        if incoming
            [selected,P,error] = nrccie_incoming_skeleton(source,A,local_tol);
        else
            [selected,P,error] = nrccie_outgoing_skeleton(source,A,local_tol);
        end
        id_seconds = id_seconds+toc(id_timer);
    end

    function L = weighted_map(source,P,selected)
        if isempty(selected)
            L = P;
            return
        end
        w = sqrt(repelem(source.wts(:),3));
        if incoming
            L = (w.*P)./w(selected).';
        else
            L = (w(selected).*P)./w.';
        end
    end
end

function bound = gram_bound(G)
% ||L||_2 <= ||L||_F = sqrt(trace(L'*L)) (or trace(L*L')).
bound = sqrt(max(0,real(trace(G))));
end

function record = node_record(children,level,nparts,ninput,nselected,tol,error,bound,G,seconds)
record = struct('children',children,'level',level,'part_count',nparts, ...
    'input_nodes',ninput,'skeleton_nodes',nselected,'requested_tol',tol, ...
    'local_error',error,'error_norm_bound',bound, ...
    'interpolation_norm_bound',gram_bound(G),'seconds',seconds);
end

function T = node_subset(S,nodes)
T = struct('npts',numel(nodes),'r',S.r(:,nodes),'n',S.n(:,nodes), ...
    'dru',S.dru(:,nodes),'drv',S.drv(:,nodes),'wts',S.wts(nodes));
end
