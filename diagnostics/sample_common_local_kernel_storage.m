% Measure actual kernel sparsity on small target samples of dominant blocks.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
[geometry,metadata] = build_leaf_geometry(6);
estimate = load(fullfile(root,'diagnostics','common_local_memory_estimate.mat'));
assert(isfield(estimate,'geometry_revision') && ...
    strcmp(estimate.geometry_revision,metadata.geometry_revision), ...
    'Stale memory estimate: run diagnostics/estimate_common_local_memory.m first.');
pattern_entries = estimate.pattern_entries; nodes = estimate.nodes; orders = estimate.orders;
actual = reshape([repmat(geometry.common.npts,1,4); ...
    cellfun(@(s)s.npts,geometry.local_tips)],1,[]);
assert(isequal(nodes(orders==6,:),actual), ...
    'Memory estimate node counts do not match current geometry.');
fprintf('order12 scalar pattern by block, billions:\n');disp(pattern_entries(:,:,end)/1e9)
zk = 2*pi/(10*69.25); eps_quad = 1e-7;
sample_count = 32; rng(83)
masters = [{geometry.common},geometry.local_tips];
samples = zeros(5,6);
capacity_samples = zeros(5,4,4); % raw nzmax, packed nzmax, raw bytes, packed bytes
for j = 1:5
    S = masters{j};
    ids = sort(randperm(S.npts,sample_count));
    T = struct('r',S.r(:,ids),'patch_id',S.patch_id(ids),'uvs_targ',S.uvs_targ(:,ids));
    B = cell(1,4); timer = tic;
    [B{1},Q] = em3d.slp.get_quad_corr_mat(S,eps_quad,zk,T);
    [B{2:4}] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk,T);
    assert(all(cellfun(@(b)all(isfinite(nonzeros(b))),B)))
    samples(j,:) = [Q.nquad,cellfun(@nnz,B),sum(cellfun(@nnz,B))/(4*Q.nquad)];
    for k = 1:4
        b = B{k}; wb = whos('b');
        transposed = b.';
        assert(nzmax(transposed)==nnz(b))
        chunks = {transposed(:,1:16),transposed(:,17:end)};
        packed = horzcat(chunks{:}).'; wp = whos('packed');
        capacity_samples(j,k,:) = [nzmax(b),nzmax(packed),wb.bytes,wp.bytes];
        assert(nzmax(packed)==nnz(packed))
        assert(wp.bytes==24*nnz(packed)+8*(size(packed,2)+1))
    end
    fprintf('master %d: pattern %d; nnz %d %d %d %d; retained %.4f; %.1f seconds\n', ...
        j,samples(j,:),toc(timer))
    fprintf('  max raw capacity/nnz %.3f; packed capacity/nnz 1\n', ...
        max(capacity_samples(j,:,1)./samples(j,2:5)))
end
geometry_revision = metadata.geometry_revision;
save(fullfile(root,'diagnostics','common_local_kernel_storage_sample.mat'), ...
    'geometry_revision', ...
    'samples','sample_count','capacity_samples')
disp('KERNEL_STORAGE_SAMPLE_PASS')
