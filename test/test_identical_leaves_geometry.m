% Geometry-only regression: four rigid copies, C4 symmetry and z reflection.
clear
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')))
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
clear pth dir
opts = struct('norder',1,'cap_mesh_spacing',0.5*69.25, ...
    'wall_profile_refinement',1);
[S,parts] = generate_leaves_surfer(.6925,opts);
base = parts.base_surfer;
assert(S.npts==4*base.npts && S.npatches==4*base.npatches)
assert(isequal([parts.node_ids{:}],1:S.npts))
assert(isequal([parts.patch_ids{:}],1:S.npatches))
assert(isequal(parts.leaf_ids,1:4))
assert(strcmp(parts.geometry_revision,leaf_geometry_revision()))
assert(all(base.wts>0) && sum(sum(base.r.*base.n,1).*base.wts(:).')>0)

scale = max(1,max(abs(S.r),[],'all'));
position_tol = 2e4*eps(scale);
unit_tol = 2e4*eps;
weight_tol = 2e4*eps(max(abs(base.wts)));
cycle = [1,2,4,3];
Q = [0,1,0;-1,0,0;0,0,1];
[reflection,upper_nodes] = leaf_midplane_permutation(base,parts.base_parts);
assert(numel(upper_nodes)==base.npts/2)
assert(isequal(reflection(reflection),1:base.npts))
fields = {'r','n','dru','drv'};
for k = 1:4
    leaf = parts.leaves{k};
    R = parts.rotation(:,:,k);
    ids = parts.node_ids{k}; patches = parts.patch_ids{k};
    assert(leaf.npts==base.npts && leaf.npatches==base.npatches)
    assert(isequal(leaf.norders,base.norders) && isequal(leaf.iptype,base.iptype))
    assert(isequal(S.norders(patches),leaf.norders) && isequal(S.iptype(patches),leaf.iptype))
    assert(isequal(S.ixyzs(patches)-ids(1)+1,leaf.ixyzs(1:end-1)))
    assert(max(abs(leaf.r-(R*base.r+parts.shift(:,k))),[],'all')<position_tol)
    assert(max(abs(leaf.wts-base.wts))<weight_tol)
    assert(max(abs(S.wts(ids)-leaf.wts))<weight_tol)
    for j = 1:numel(fields)
        field = fields{j};
        assert(max(abs(S.(field)(:,ids)-leaf.(field)),[],'all')<position_tol)
        if j>1
            assert(max(abs(leaf.(field)-R*base.(field)),[],'all')<unit_tol)
        end
    end
    next = parts.leaves{cycle(mod(find(cycle==k),4)+1)};
    assert(max(abs(Q*leaf.r-next.r),[],'all')<position_tol)
    for j = 2:numel(fields)
        field = fields{j};
        assert(max(abs(Q*leaf.(field)-next.(field)),[],'all')<unit_tol)
    end
    assert(max(abs(leaf.r(:,reflection)-[1;1;-1].*leaf.r),[],'all')<position_tol)
    assert(max(abs(leaf.n(:,reflection)-[1;1;-1].*leaf.n),[],'all')<unit_tol)
end
fprintf('IDENTICAL_LEAVES_GEOMETRY_PASS: %d nodes / leaf, %d total.\n',base.npts,S.npts)
