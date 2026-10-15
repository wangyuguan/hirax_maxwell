function [F,info] = nrccie_srskelf(S,leaves,C4,opts,max_nodes,tol)
%
%  nrccie_srskelf
%    Factor the corrected four-leaf NRCCIE operator using RS-S.
%
%  Input arguments:
%    * S: merged surfer; leaves contains its four components
%    * C4: compact quadrature corrections from fourfold_quad_corr_mats
%    * opts: zk, alpha, jump; optional tree_axes ([1,2,3]) and verb (0)
%    * max_nodes: physical-node occupancy of each RS-S box
%    * tol: local ID tolerance; true residuals are checked by the solver
%
%  Output arguments:
%    * F: asymmetric strong-skel factorization
%    * info: factor sizes, timings, and tree parameters
%

    validateattributes(max_nodes,{'numeric'}, ...
        {'scalar','real','finite','integer','positive'})
    validateattributes(tol,{'numeric'}, ...
        {'scalar','real','finite','positive','<',1})
    assert(iscell(leaves) && numel(leaves)==4 && ...
        S.npts==sum(cellfun(@(leaf)leaf.npts,leaves)), ...
        'Expected S to be the merge of four leaves.')
    axes = [1,2,3];
    if isfield(opts,'tree_axes')
      axes = opts.tree_axes;
    end
    validateattributes(axes,{'numeric'}, ...
        {'vector','nonempty','real','finite','integer','>=',1,'<=',3})
    assert(numel(unique(axes))==numel(axes),'tree_axes must be distinct.')
    coordinates = repelem(S.cms(axes,S.patch_id),1,3);
    assert(size(coordinates,2)==3*S.npts && all(isfinite(coordinates(:))), ...
        'Expected finite patch centers for every surface node.')
    occupancy = 3*max_nodes;
    [~,~,groups] = unique(coordinates.','rows');
    largest_group = max(accumarray(groups,1));
    assert(largest_group<=occupancy, ...
        ['Tree occupancy %d cannot split %d coincident scalar DOFs; ' ...
        'increase max_nodes or add a tree axis.'], ...
        occupancy,largest_group)

    assert(exist('srskelf_asym','file')==2 && exist('id','file')==2, ...
        'Strong-skel and FLAM must be available; run startup first.')
    factor_options = struct('symm','n','verb',0);
    if isfield(opts,'verb')
      factor_options.verb = opts.verb;
    end
    matrix_block = @(rows,cols)nrccie_fourfold_system_block( ...
        leaves,C4,opts,rows,cols);
    timer = tic;
    F = srskelf_asym(matrix_block,coordinates,occupancy,tol,[],factor_options);
    seconds = toc(timer);

    eliminated = reshape(arrayfun(@(box)numel(box.rd),F.factors),1,[]);
    prefix = [0,cumsum(eliminated)];
    remaining_by_level = F.N-prefix(F.lvp+1);
    root_slots = F.lvp(end-1)+1:F.lvp(end);
    storage = whos('F');
    info = struct('method','direct RS-S of complete four-leaf A', ...
        'local_tolerance',tol,'occupancy',occupancy,'tree_axes',axes, ...
        'levels',F.nlvl,'remaining_by_level',remaining_by_level, ...
        'root_dofs',sum(eliminated(root_slots)), ...
        'factor_bytes',storage.bytes,'seconds',seconds, ...
        'randomized',false,'proxy',false,'full_matrix_materialized',false);
    if factor_options.verb
      fprintf('Direct A RS-S: root %d/%d DOFs, %.3f GiB, %.2f s.\n', ...
          info.root_dofs,F.N,info.factor_bytes/2^30,seconds)
    end
end
