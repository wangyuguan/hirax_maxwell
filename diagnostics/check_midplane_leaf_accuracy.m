% Test mirror covariance on true order-6 leaf patches, without a full matrix.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
[geometry,metadata] = build_leaf_geometry(6);
rng(26); zk = 2*pi/(10*69.25); tolerances = [1e-9,1e-11];
source_info.r = [0;0;10*69.25]; source_info.edips = -69.25^3*[1;1i;0];
matrix_errors = zeros(5,2,4); smooth_errors = matrix_errors;
signs = [1,1,1,-1];
for j = 1:5
    if j==1
        S = geometry.common; g = metadata.common_groups; leaf = 1;
    else
        leaf = j-1; S = geometry.local_tips{leaf}; g = metadata.local_groups{leaf};
    end
    S = affine_transf(S,geometry.rotation(:,:,leaf),geometry.shift(:,leaf));
    [perm,upper] = leaf_midplane_permutation(S,g);
    % Stratify the samples over core, collar, and wall rather than area alone.
    ids = [];
    for name = {'top_core','top_collar','upper_wall'}
        patches = g.(name{1}); group = find(ismember(S.patch_id,patches));
        ids = [ids,reshape(group(randperm(numel(group),min(24,numel(group)))),1,[])]; %#ok<AGROW>
    end
    ids = sort(ids); assert(all(ismember(ids,upper)))
    [ei,hi] = em3d.incoming_sources(zk,source_info,S,'electric dipole');
    density = [cross(S.n,hi,1);sum(S.n.*ei,1)].';
    for it = 1:numel(tolerances)
        B = nrccie_quad_corr_block(S,tolerances(it),zk,[],72,[ids,perm(ids)]);
        n = numel(ids);
        for k = 1:4
            direct = B{k}(n+1:end,:);
            reflected = signs(k)*B{k}(1:n,perm);
            matrix_errors(j,it,k) = norm(direct-reflected,'fro')/norm(direct,'fro');
            v = direct*density;
            smooth_errors(j,it,k) = norm(v-reflected*density,'fro')/norm(v,'fro');
        end
        fprintf('master %d eps %.0e: matrix %.3e %.3e %.3e %.3e; incident-trace action %.3e %.3e %.3e %.3e\n', ...
            j,tolerances(it),matrix_errors(j,it,:),smooth_errors(j,it,:))
    end
end
geometry_revision = metadata.geometry_revision;
save(fullfile(root,'diagnostics','midplane_leaf_accuracy.mat'), ...
    'geometry_revision', ...
    'matrix_errors','smooth_errors','tolerances','zk')
assert(all(isfinite(matrix_errors),'all') && all(isfinite(smooth_errors),'all'))
assert(max(smooth_errors(:,end,:),[],'all')<5e-9, ...
    'Mirror consistency on incident traces is insufficient for a 1e-8 solve.')
disp('MIDPLANE_LEAF_ACCURACY_PASS')
