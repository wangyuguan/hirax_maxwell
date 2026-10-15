% Compare original and advanced straight partitions on identical outlines.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')))
cases = [0 2; 1 2; 1 4]; % promote identical connectors, interface spacing
results = zeros(size(cases,1),8);
for j = 1:size(cases,1)
    opts = struct('extend_common_outline',logical(cases(j,1)), ...
        'interface_spacing',cases(j,2));
    [~,~,S,p] = build_leaf_geometry(4,opts);
    if j==1
        reference = p.outlines;
    else
        for leaf = 1:4
            a = reference{leaf}.position_coefficients;
            b = p.outlines{leaf}.position_coefficients;
            index = find(arrayfun(@(k)isequal(a(:,:,1),b(:,:,k)),1:size(b,3)),1);
            assert(~isempty(index) && isequal(a,circshift(b,1-index,3)), ...
                'Repartitioning changed the physical outline.');
        end
    end
    common = p.local_common_surfer.npts;
    local = sum(cellfun(@(s)s.npts,p.local_unique_surfers));
    ca = sum(p.local_common_surfer.wts);
    la = sum(cellfun(@(s)sum(s.wts),p.local_unique_surfers));
    results(j,:) = [cases(j,:),S.npts,common,local,common+local, ...
        4*ca/(4*ca+la),size(p.interface_nodes,2)];
    fprintf('promote=%d spacing=%.1f total=%d common=%d local=%d stored=%d area=%.6f interface=%d\n',results(j,:));
end
geometry_revision = leaf_geometry_revision();
save(fullfile(root,'diagnostics','common_partition_comparison.mat'), ...
    'results','cases','geometry_revision');
