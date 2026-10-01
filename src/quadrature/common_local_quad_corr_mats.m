function [C,info] = common_local_quad_corr_mats(common,local,eps_quad,zk,batch_size,reflection)
% Four world-frame common/local pairs, ordered [U1,D1,U2,D2,U3,D3,U4,D4].
% Only U-U blocks are rotationally reused. All blocks involving D are exact.
% Optional REFLECTION is a cell of eight midplane node permutations.
if nargin<5, batch_size = 2000; end
if nargin<6, reflection = []; end
common_reflection = [];
if ~isempty(reflection), common_reflection = reflection{1}; end
[C.common,common_info] = fourfold_quad_corr_mats( ...
    common{1},common{2},common{3},common{4},eps_quad,zk,batch_size,common_reflection);
pieces = reshape([common(:).';local(:).'],1,[]);
counts = cellfun(@(s)s.npts,pieces);
offset = [0,cumsum(counts)];
C.nodes = cell(1,8);
for j = 1:8
    C.nodes{j} = offset(j)+(1:counts(j));
end
C.common_nodes = [C.nodes{1:2:8}];
C.npts = offset(end);
C.blocks = cell(8,8);
C.reflection = reflection;
local_nonzeros = 0;
for target = 1:8
    target_ids = [];
    if ~isempty(reflection)
        target_ids = find((1:counts(target))<reflection{target});
    end
    for source = 1:8
        if mod(target,2)==1 && mod(source,2)==1, continue; end
        T = pieces{target};
        if target==source, T = []; end
        B = nrccie_quad_corr_block(pieces{source},eps_quad,zk,T,batch_size,target_ids);
        n = sum(cellfun(@nnz,B));
        if n>0, C.blocks{target,source} = B; end
        local_nonzeros = local_nonzeros+n;
    end
end
info.common = common_info;
info.local_nonzeros = local_nonzeros;
info.stored_nonzeros = common_info.stored_nonzeros+local_nonzeros;
info.assembled_nonzeros = common_info.assembled_nonzeros+(1+~isempty(reflection))*local_nonzeros;
info.midplane_reused = ~isempty(reflection);
w = whos('C');
info.stored_bytes = w.bytes;
end
