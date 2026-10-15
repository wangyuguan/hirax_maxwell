% Estimate node counts and correction storage without evaluating quadrature.
% Counts are exact for the current generated patch partition. Near patterns are sampled
% on its order-6 polynomial geometry evaluated at each requested node set.
% Sparse bytes assume 16-byte complex values + 8-byte row indices per entry;
% exact zero entries could reduce storage. Column pointers are not included.
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
[geometry,metadata] = build_leaf_geometry(6);
orders = [4,6,8,10,12];
sample_count = 8192;
gmres_restart = 100;
rng(29)
masters = [{geometry.common},geometry.local_tips];
patch_counts = zeros(5,2);
for j = 1:5
    patch_counts(j,:) = [nnz(masters{j}.iptype==1),nnz(masters{j}.iptype==11)];
end
expanded_patch_counts = 4*patch_counts(1,:)+sum(patch_counts(2:end,:),1);
fprintf('Master patch counts [triangles, quads]:\n'); disp(patch_counts)
fprintf('Expanded patch counts: %d triangles, %d quads\n',expanded_patch_counts)
fprintf('Actual order6 master bytes: %.3f GiB\n',whos_bytes('geometry')/2^30)

sources = cell(1,8);
for j = 1:4
    for part = 1:2
        m = 1+(part-1)*j;
        sources{2*j-2+part} = struct('cms',geometry.rotation(:,:,j)*masters{m}.cms+ ...
            geometry.shift(:,j),'rads',masters{m}.rads,'iptype',masters{m}.iptype);
    end
end
nodes = zeros(numel(orders),8);
pattern_entries = zeros(8,8,numel(orders));
pattern_se = zeros(8,8,numel(orders));
memory_gib = zeros(numel(orders),4);
for k = 1:numel(orders)
    p = orders(k); ntri = (p+1)*(p+2)/2; nquad = (p+1)^2;
    rfac = 1.25; if p<=6, rfac = 2; end
    vm = {koorn.coefs2vals(6,koorn.rv_nodes(p)), ...
        polytens.lege.coefs2vals(6,polytens.lege.nodes(p))};
    for target = 1:8
        leaf = ceil(target/2); m = 1+mod(target+1,2)*leaf;
        S = masters{m};
        sizes = ntri*(S.iptype==1)+nquad*(S.iptype==11);
        offsets = [0;cumsum(sizes(:))]; nodes(k,target) = offsets(end);
        ids = sort(randperm(offsets(end),min(sample_count,offsets(end))));
        % histcounts edges are zero-based while MATLAB node IDs are one-based.
        [~,~,patch] = histcounts(ids-1,[offsets;Inf]);
        points = zeros(3,numel(ids));
        for j = unique(patch)
            chosen = find(patch==j); local_id = ids(chosen)-offsets(j);
            points(:,chosen) = S.srccoefs{j}(1:3,:)*vm{1+(S.iptype(j)==11)}(local_id,:).';
        end
        points = geometry.rotation(:,:,leaf)*points+geometry.shift(:,leaf);
        for source = 1:8
            src = sources{source};
            [row,col] = near_indices(src.cms,src.rads*rfac,points);
            srcsizes = ntri*(src.iptype==1)+nquad*(src.iptype==11);
            prefix = [0;cumsum(srcsizes(col))];
            per_target = prefix(row(2:end))-prefix(row(1:end-1));
            pattern_entries(target,source,k) = mean(per_target)*nodes(k,target);
            pattern_se(target,source,k) = std(per_target)/sqrt(numel(ids))*nodes(k,target);
        end
    end
    a = pattern_entries(:,:,k); se = pattern_se(:,:,k);
    uu = a(1,1:2:8); a(1:2:8,1:2:8) = 0;
    e = se(1,1:2:8); se(1:2:8,1:2:8) = 0;
    stored_pattern = sum(uu)+sum(a,'all');
    largest = max([uu(:);a(:)]);
    running = 0; build_peak = 0;
    ordered = [uu([1,2,4,3]),reshape(a.',1,[])];
    for n = ordered
        % Existing chunks plus two temporary copies of one scalar kernel.
        build_peak = max(build_peak,running+144*n);
        running = running+96*n;
    end
    memory_gib(k,:) = [96*stored_pattern,48*gmres_restart*sum(nodes(k,:)), ...
        build_peak,96*sum(pattern_entries(:,:,k),'all')]/2^30;
    fprintf(['order %2d: N=%d; unknowns=%d; stored scalar pattern=%.6g +/- %.2g; ' ...
        'four complex sparse upper estimate=%.2f GiB; largest block pattern=%.6g; ' ...
        'GMRES%d basis=%.2f GiB\n'],p,sum(nodes(k,:)),3*sum(nodes(k,:)), ...
        stored_pattern,2*sqrt(sum(e.^2)+sum(se.^2,'all')), ...
        memory_gib(k,1),largest,gmres_restart,memory_gib(k,2))
end
disp(array2table([orders(:),sum(nodes,2),memory_gib], ...
    'VariableNames',{'order','nodes','stored_GiB','GMRES_basis_GiB', ...
    'batch_concat_peak_GiB','expanded_GiB'}))
midplane_memory_gib = memory_gib;
midplane_memory_gib(:,[1,3]) = midplane_memory_gib(:,[1,3])/2;
disp('With midplane reflection (basis and full expanded operator unchanged):')
disp(array2table([orders(:),sum(nodes,2),midplane_memory_gib], ...
    'VariableNames',{'order','nodes','stored_GiB','GMRES_basis_GiB', ...
    'batch_concat_peak_GiB','expanded_GiB'}))
geometry_revision = metadata.geometry_revision;
save(fullfile(root,'diagnostics','common_local_memory_estimate.mat'), ...
    'geometry_revision', ...
    'orders','nodes','patch_counts','expanded_patch_counts','sample_count', ...
    'pattern_entries','pattern_se','memory_gib','midplane_memory_gib','gmres_restart')
disp('MEMORY_ESTIMATE_PASS')

function [row,col] = near_indices(cms,radii,targs)
np = size(cms,2); nt = size(targs,2);
id = 'findnearmem(i double[xx], i int64_t[x], i double[x], i int64_t[x], i double[xx], i int64_t[x], io int64_t[x])';
n = fmm3dbie_routs(id,cms,np,radii,3,targs,nt,0,3,np,1,np,1,3,nt,1,1);
id = 'findnear(i double[xx], i int64_t[x], i double[x], i int64_t[x], i double[xx], i int64_t[x], io int64_t[x], io int64_t[x])';
[row,col] = fmm3dbie_routs(id,cms,np,radii,3,targs,nt,zeros(nt+1,1),zeros(n,1), ...
    3,np,1,np,1,3,nt,1,nt+1,n);
end

function n = whos_bytes(name)
w = evalin('caller',sprintf('whos(''%s'')',name)); n = w.bytes;
end
