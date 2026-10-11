function A = nrccie_fourfold_system_block(leaves,C4,opts,rows,cols)
%NRCCIE_FOURFOLD_SYSTEM_BLOCK Selected entries of the corrected four-leaf A.
% Scalar DOFs are in merged-leaf order and interleaved [ju,jv,rho] at
% every node. C4 symmetry reduces every leaf pair to one of the four
% interactions with leaf 1. The complete matrix is never assembled.

assert(iscell(leaves) && numel(leaves)==4, ...
    'nrccie_fourfold_system_block:Leaves','Expected four leaves.')
assert(isstruct(C4) && isequal(C4.cycle_to_merge,[1,2,4,3]), ...
    'nrccie_fourfold_system_block:Cycle', ...
    'Expected the C4 rotation cycle [1,2,4,3].')
assert(~isfield(C4,'reflection') || isempty(C4.reflection), ...
    'nrccie_fourfold_system_block:ReflectionUnsupported', ...
    'Reflection-compressed corrections are not supported here.')

npts = C4.component_npts;
assert(all(cellfun(@(leaf)leaf.npts==npts,leaves)), ...
    'nrccie_fourfold_system_block:Size', ...
    'Leaves and C4 corrections must have matching node counts.')
n = 3*npts;
N = 4*n;
rows = checked_indices(rows,N);
cols = checked_indices(cols,N);
if numel(rows)>=N && numel(cols)>=N && ...
        numel(unique(rows))==N && numel(unique(cols))==N
    error('nrccie_fourfold_system_block:FullMatrix', ...
        'Refusing to assemble the complete four-leaf matrix A.')
end

A = complex(zeros(numel(rows),numel(cols)));
if isempty(rows) || isempty(cols), return; end

row_leaf = floor((rows-1)/n)+1;
col_leaf = floor((cols-1)/n)+1;
local_rows = mod(rows-1,n)+1;
local_cols = mod(cols-1,n)+1;
cycle = C4.cycle_to_merge;

for target_cycle = 1:4
    ri = find(row_leaf==cycle(target_cycle));
    if isempty(ri), continue; end
    for source_cycle = 1:4
        ci = find(col_leaf==cycle(source_cycle));
        if isempty(ci), continue; end
        relative = mod(source_cycle-target_cycle,4)+1;
        C = {C4.slp{relative},C4.gx{relative}, ...
            C4.gy{relative},C4.gz{relative}};
        target = leaves{1};
        if relative==1, target = []; end
        A(ri,ci) = nrccie_system_block(leaves{cycle(relative)},target, ...
            C,opts,local_rows(ri),local_cols(ci));
    end
end
end


function indices = checked_indices(indices,N)
assert(isnumeric(indices) && (isvector(indices) || isempty(indices)) && ...
    isreal(indices) && all(isfinite(indices(:))) && ...
    all(indices(:)==fix(indices(:))) && ...
    all(indices(:)>=1 & indices(:)<=N), ...
    'nrccie_fourfold_system_block:Indices', ...
    'Expected scalar DOF indices in 1:N.')
indices = reshape(indices,1,[]);
end
