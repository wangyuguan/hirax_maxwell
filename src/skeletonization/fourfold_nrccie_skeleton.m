function [J,blocks,Pmerge,info] = fourfold_nrccie_skeleton(leaves,samplefun,tol,max_patches,direction)
%FOURFOLD_NRCCIE_SKELETON One-leaf multilevel ID reused on four rotated copies.
%
% Input arguments:
%   leaves      - four surfers in merge order [1 2 3 4], with clockwise
%                 rotation cycle [1 2 4 3] about the z axis through the origin.
%   samplefun   - @(T) actual samples: 3*T.npts rows for incoming or columns
%                 for outgoing. The fixed sample set must be C4-covariant:
%                 rotating a leaf changes samples only by a unitary map on
%                 the sample dimension (e.g. a weighted sample permutation).
%   tol         - relative weighted Frobenius sample-error tolerance.
%   max_patches - maximum patches per finest-level connected region.
%   direction   - 'incoming' or 'outgoing'.
%
% Output arguments:
%   J, blocks, Pmerge - same maps as MULTILEVEL_NRCCIE_SKELETON, in full
%                       four-leaf merge order. Pmerge is sparse block diagonal.
%   info - representative-leaf hierarchy/timings, reuse and storage metadata.
%
% IDs are computed only on leaves{1}; no ID merges different leaves. Matched
% rotated tangent frames preserve the scalar coordinates [ju,jv,rho], so all
% four leaves use the same interpolation coefficients. The inherited error
% bound applies to the four-leaf samples under the covariance assumption;
% this routine checks geometry, not sample covariance or unseen fields.

assert(iscell(leaves) && numel(leaves)==4, ...
    'Supply four surfaces in merge order [1 2 3 4].')
assert(any(strcmp(direction,{'incoming','outgoing'})))
timer = tic;
cycle = [1,2,4,3];
Q = [0,1,0;-1,0,0;0,0,1];
S = leaves{1};
npts = S.npts;
fields = {'r','n','dru','drv','wts'};
symmetry_errors = zeros(4,numel(fields));
for leaf = 1:4
    T = leaves{leaf};
    assert(T.npts==npts && T.npatches==S.npatches && ...
        isequal(T.norders,S.norders) && isequal(T.iptype,S.iptype) && ...
        isequal(T.ixyzs,S.ixyzs) && isequal(T.patch_id,S.patch_id), ...
        'Rotational reuse requires identical node and patch ordering.')
    assert(all(T.wts>0),'Skeletonization requires positive source weights.')
    power = find(cycle==leaf)-1;
    rotation = Q^power;
    for k = 1:numel(fields)
        name = fields{k};
        expected = S.(name);
        if ~strcmp(name,'wts'), expected = rotation*expected; end
        actual = T.(name);
        scale = max(1,max(abs(expected),[],'all'));
        symmetry_errors(leaf,k) = max(abs(expected-actual),[],'all');
        assert(symmetry_errors(leaf,k)<=1e3*eps(scale), ...
            'Rotational reuse requires C4-matched positions, frames and weights.')
    end
end
geometry_seconds = toc(timer);

[local_J,local_blocks,local_Pmerge,local_info] = ...
    multilevel_nrccie_skeleton(S,samplefun,tol,max_patches,direction);
lift_timer = tic;
nblocks = numel(local_blocks);
ncandidate = 3*numel(local_info.candidate_nodes);
J = reshape(local_J(:)+3*npts*(0:3),1,[]);
blocks = repmat(local_blocks,4,1);
for leaf = 1:4
    for p = 1:nblocks
        block = (leaf-1)*nblocks+p;
        blocks(block).nodes = local_blocks(p).nodes+(leaf-1)*npts;
        blocks(block).skeleton_nodes = ...
            local_blocks(p).skeleton_nodes+(leaf-1)*npts;
        blocks(block).candidate_dofs = ...
            local_blocks(p).candidate_dofs+(leaf-1)*ncandidate;
    end
end
Pmerge = kron(speye(4),sparse(local_Pmerge));

info = local_info;
info.method = 'multilevel-c4-reuse';
info.reused = true;
info.representative_leaf = 1;
info.representative_J = local_J;
info.representative_info = local_info;
info.cycle = cycle;
info.nodes_per_leaf = npts;
info.skeleton_nodes_per_leaf = (numel(local_J)/3)*ones(1,4);
info.block_count_per_leaf = nblocks;
info.candidate_dofs_per_leaf = ncandidate;
info.symmetry_fields = fields;
info.symmetry_errors = symmetry_errors;
info.tree_scope = 'representative leaf only';
info.error_bound_scope = 'derived for C4-covariant supplied samples';
info.reuse_sample_error_validated = false;
info.part_id = repmat(local_info.part_id,4,1)+ ...
    repelem((0:3).'*nblocks,S.npatches);
info.adjacency = kron(speye(4),local_info.adjacency);
info.part_patch_counts = repmat(local_info.part_patch_counts,4,1);
info.part_node_counts = repmat(local_info.part_node_counts,4,1);
info.local_skeleton_counts = repmat(local_info.local_skeleton_counts,4,1);
info.local_errors = repmat(local_info.local_errors,4,1);
info.candidate_nodes = reshape(local_info.candidate_nodes(:)+npts*(0:3),1,[]);
info.full_sample_mib = 4*local_info.full_sample_mib;
storage = whos('blocks','Pmerge');
info.interpolation_mib = sum([storage.bytes])/2^20;
info.geometry_seconds = geometry_seconds;
info.reuse_seconds = toc(lift_timer);
info.seconds = toc(timer);
end
