function [nodes,U,info] = maxwell_multilevel_row_id(P,rowfun,tol,max_leaf_points)
% MAXWELL_MULTILEVEL_ROW_ID Deterministic joint ID on sphere point blocks.
%
%  Input arguments:
%    * P: surface geometry with P.r of size 3 by P.npts
%    * rowfun: callback F = rowfun(point_ids), returning four consecutive
%        rows per requested point, in exactly the requested order. For the
%        joint incoming/outgoing ID, return
%          [B(:,dofs)'/normB, C(dofs,:)/normC],
%        where dofs = reshape(4*(point_ids-1)+(1:4).',1,[]). The prime is
%        the conjugate transpose. Use fixed normalizations for all calls.
%    * tol: requested relative Frobenius tolerance (default 1e-6)
%    * max_leaf_points: maximum points in a spatial leaf (default 32)
%
%  Output arguments:
%    * nodes: root skeleton point indices, in ascending order
%    * U: full interpolation matrix for this small-problem implementation:
%          F approximately equals U*F(info.dofs,:).
%        All four components at a selected point are retained, and
%        U(info.dofs,:) is exactly the identity. Thus B approximately
%        equals B(:,info.dofs)*U', and C approximately equals
%        U*C(info.dofs,:).
%    * info: tree and local transfers, root DOFs, residuals, and statistics
%
%  The spatial tree bisects the longest coordinate extent at its median.
%  Leaves compress their point blocks. A parent only recompresses the union
%  of its children's skeletons. No randomized sketch or preassembled global
%  signature is required. If compression is ineffective, a candidate block
%  can be as large as the full point set at the root.
%
%  Each local ID uses tol/nlevels. This allocation is a heuristic because
%  interpolation can amplify child errors. info.relative_error_bound is a
%  conservative bound including that amplification. A final streaming pass
%  over the leaves measures info.relative_residual against the actual joint
%  signature; this adds callback evaluations but does not store global F.
%  The requested tolerance is not an unconditional global guarantee.
%
%  info.tree(k).transfer interpolates candidate rows from the selected
%  skeleton rows at box k. Together with the children and point indices, these
%  transfers allow a future implementation to apply U implicitly. Here the
%  full U is returned to simplify validation of small scattering experiments.

    if nargin < 3 || isempty(tol), tol = 1e-6; end
    if nargin < 4 || isempty(max_leaf_points), max_leaf_points = 32; end
    assert(isa(rowfun,'function_handle'),'rowfun must be a function handle.');
    assert(isscalar(tol) && isreal(tol) && isfinite(tol) && tol>0 && tol<1, ...
        'tol must lie between zero and one.');
    assert(isscalar(max_leaf_points) && isfinite(max_leaf_points) && ...
        max_leaf_points>=1 && max_leaf_points==fix(max_leaf_points), ...
        'max_leaf_points must be a positive integer.');
    npts = P.npts;
    xyz = P.r;
    assert(isscalar(npts) && npts>=1 && npts==fix(npts), ...
        'P.npts must be a positive integer.');
    assert(isequal(size(xyz),[3,npts]) && isreal(xyz) && ...
        all(isfinite(xyz(:))),'P.r must be a finite real 3 by P.npts array.');

    empty_box = struct('point_ids',[],'children',[],'depth',0, ...
        'candidate_nodes',[],'skeleton_nodes',[],'transfer',[], ...
        'relative_residual',0,'absolute_residual',0,'error_bound',0, ...
        'interpolation_norm_bound',0,'rank',0,'fallback',false);
    tree = repmat(empty_box,0,1);
    root = build_tree(1:npts,1);
    nlevels = max([tree.depth]);
    local_tol = tol/nlevels;
    callback_calls = 0;
    max_block_points = 0;
    signature_columns = [];
    [nodes,U,error_bound] = compress_box(root);

    % Validate by streaming original leaf rows. Only the root skeleton
    % signature and one leaf block need to coexist during this pass.
    if isempty(nodes)
        root_rows = zeros(0,signature_columns,'like',U);
    else
        root_rows = get_rows(nodes);
    end
    leaf_boxes = find(arrayfun(@(box) isempty(box.children),tree));
    residual_squared = 0;
    scale_squared = 0;
    for k = leaf_boxes(:).'
        check_ids = tree(k).point_ids;
        check_rows = get_rows(check_ids);
        check_residual = check_rows-U(point_dofs(check_ids),:)*root_rows;
        residual_squared = residual_squared+norm(check_residual,'fro')^2;
        scale_squared = scale_squared+norm(check_rows,'fro')^2;
    end
    scale = sqrt(scale_squared);
    if scale==0
        relative_residual = 0;
        relative_error_bound = 0;
    else
        relative_residual = sqrt(residual_squared)/scale;
        relative_error_bound = error_bound/scale;
    end
    if isempty(U), maxabsU = 0; else, maxabsU = max(abs(U(:))); end
    info = struct('dofs',point_dofs(nodes),'npoints',numel(nodes), ...
        'rank',tree(root).rank,'nlevels',nlevels,'nboxes',numel(tree), ...
        'nleaves',numel(leaf_boxes),'local_tolerance',local_tol, ...
        'relative_residual',relative_residual, ...
        'relative_error_bound',relative_error_bound, ...
        'max_local_relative_residual',max([tree.relative_residual]), ...
        'maxabsU',maxabsU,'callback_calls',callback_calls, ...
        'max_block_points',max_block_points, ...
        'fallback_boxes',find([tree.fallback]),'root',root,'tree',tree);

    function box_id = build_tree(ids,depth)
        box_id = numel(tree)+1;
        tree(box_id) = empty_box;
        tree(box_id).point_ids = ids;
        tree(box_id).depth = depth;
        if numel(ids)<=max_leaf_points, return; end
        extent = max(xyz(:,ids),[],2)-min(xyz(:,ids),[],2);
        [~,axis_id] = max(extent);
        % The point index breaks ties deterministically, including the case
        % of coincident coordinates. Store every regional index list sorted.
        ordered = sortrows([xyz(axis_id,ids).',ids.'],[1,2]);
        midpoint = floor(numel(ids)/2);
        left_ids = sort(ordered(1:midpoint,2)).';
        right_ids = sort(ordered(midpoint+1:end,2)).';
        left = build_tree(left_ids,depth+1);
        right = build_tree(right_ids,depth+1);
        tree(box_id).children = [left,right];
    end

    function rows = get_rows(ids)
        rows = rowfun(ids);
        assert(isnumeric(rows) && ismatrix(rows) && ...
            size(rows,1)==4*numel(ids), ...
            'rowfun must return four rows per requested point.');
        assert(all(isfinite(rows(:))),'rowfun returned nonfinite values.');
        if isempty(signature_columns)
            signature_columns = size(rows,2);
        else
            assert(size(rows,2)==signature_columns, ...
                'rowfun must return the same columns on every call.');
        end
        callback_calls = callback_calls+1;
        max_block_points = max(max_block_points,numel(ids));
    end

    function [selected,interpolation,bound] = compress_box(box_id)
        ids = tree(box_id).point_ids;
        children = tree(box_id).children;
        child_nodes = cell(1,numel(children));
        child_U = cell(1,numel(children));
        child_bounds = zeros(1,numel(children));
        if isempty(children)
            candidates = ids;
            child_gain = 1;
        else
            child_gain = 0;
            for j = 1:numel(children)
                [child_nodes{j},child_U{j},child_bounds(j)] = ...
                    compress_box(children(j));
                child_gain = max(child_gain, ...
                    tree(children(j)).interpolation_norm_bound);
            end
            candidates = sort([child_nodes{:}]);
        end

        if isempty(candidates)
            % Both children represented an identically zero signature.
            selected = [];
            interpolation = zeros(4*numel(ids),0,'like',child_U{1});
            transfer = zeros(0,0,'like',interpolation);
            local_relative = 0;
            local_absolute = 0;
            numerical_rank = 0;
            fallback = false;
        else
            rows = get_rows(candidates);
            [local_nodes,transfer,local_info] = ...
                maxwell_shared_row_id(rows,local_tol);
            [local_nodes,permutation] = sort(local_nodes);
            transfer = transfer(:,point_dofs(permutation));
            selected = candidates(local_nodes);
            local_relative = local_info.relative_residual;
            numerical_rank = local_info.rank;
            fallback = local_relative>local_tol*(1+32*eps(class(rows)));
            if fallback || numel(selected)==numel(candidates)
                % A rank-deficient least-squares fit must not spoil an
                % otherwise exact no-compression case.
                selected = candidates;
                transfer = eye(4*numel(candidates),'like',rows);
                local_relative = 0;
            end
            local_absolute = local_relative*norm(rows,'fro');
            if isempty(children)
                interpolation = transfer;
            else
                interpolation = zeros(4*numel(ids),4*numel(selected), ...
                    'like',transfer);
                for j = 1:numel(children)
                    [~,row_locations] = ismember(tree(children(j)).point_ids,ids);
                    [~,candidate_locations] = ismember(child_nodes{j},candidates);
                    interpolation(point_dofs(row_locations),:) = ...
                        child_U{j}*transfer(point_dofs(candidate_locations),:);
                end
            end
            [~,selected_locations] = ismember(selected,ids);
            interpolation(point_dofs(selected_locations),:) = ...
                eye(4*numel(selected),'like',interpolation);
        end

        % For disjoint child row/column blocks, the composition operator's
        % 2-norm is bounded by the maximum of the child bounds. The standard
        % sqrt(||U||_1*||U||_inf) estimate avoids a costly spectral norm.
        bound = norm(child_bounds)+child_gain*local_absolute;
        if isempty(interpolation)
            norm_bound = 0;
        else
            norm_bound = sqrt(norm(interpolation,1)*norm(interpolation,inf));
        end
        tree(box_id).candidate_nodes = candidates;
        tree(box_id).skeleton_nodes = selected;
        tree(box_id).transfer = transfer;
        tree(box_id).relative_residual = local_relative;
        tree(box_id).absolute_residual = local_absolute;
        tree(box_id).error_bound = bound;
        tree(box_id).interpolation_norm_bound = norm_bound;
        tree(box_id).rank = numerical_rank;
        tree(box_id).fallback = fallback;
    end
end

function dofs = point_dofs(nodes)
    dofs = reshape(4*(nodes(:).'-1)+(1:4).',1,[]);
end
