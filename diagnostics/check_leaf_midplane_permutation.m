% Check all five current master pieces and their world-frame placements.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
[geometry,metadata] = build_leaf_geometry(6);
pieces = [{geometry.common},geometry.local_tips];
groups = [{metadata.common_groups},metadata.local_groups];
for j = 1:5
    S = pieces{j};
    [permutation,upper] = leaf_midplane_permutation(S,groups{j});
    position_error = max(abs(S.r(:,permutation)-[1;1;-1].*S.r),[],'all');
    normal_error = max(abs(S.n(:,permutation)-[1;1;-1].*S.n),[],'all');
    weights = S.wts(:).';
    weight_error = max(abs(weights(permutation)-weights));
    assert(all(S.r(3,upper)>0) && all(S.r(3,permutation(upper))<0))
    for leaf = 1:4
        world = affine_transf(S,geometry.rotation(:,:,leaf),geometry.shift(:,leaf));
        [world_permutation,world_upper] = leaf_midplane_permutation(world,groups{j});
        assert(isequal(permutation,world_permutation) && isequal(upper,world_upper))
    end
    fprintf('Piece %d: %d nodes; reflection errors r %.3e, n %.3e, w %.3e\n', ...
        j,S.npts,position_error,normal_error,weight_error)
end
disp('MIDPLANE_PERMUTATION_PASS')
