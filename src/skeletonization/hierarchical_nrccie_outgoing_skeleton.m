function [Jout,blocks,Pmerge,info] = hierarchical_nrccie_outgoing_skeleton(S,T,zk,tol,max_patches)
%HIERARCHICAL_NRCCIE_OUTGOING_SKELETON Connected source parts and a shared proxy.
% T is a separated, closed enclosing proxy with positions r and area weights
% wts. Each part and the union of its candidates use standard QR column ID,
% completed to whole source nodes. Original source weights are preserved.
% Jout contains global interleaved [ju;jv;rho] DOFs of the final nodes.
% For a 3-by-S.npts density, apply the factored interpolation using
%   qout = zeros(numel(Jout),1);
%   for i = 1:numel(blocks)
%       density_i = reshape(density(:,blocks(i).nodes),[],1);
%       qout = qout + Pmerge(:,blocks(i).candidate_dofs) ...
%           * (blocks(i).P*density_i);
%   end
% Proxy matrices are built for one part or the candidate union at a time;
% the global interpolation is kept as local factors and Pmerge.
% info.proxy_error_bound bounds the source/target area-weighted Frobenius
% proxy error; independently validate fields at new exterior targets.

[part_id,adjacency] = partition_surface_patches(S,max_patches);
nparts = max(part_id);
node_part = part_id(S.patch_id);
proxy_w = sqrt(repelem(T.wts(:),6));
blocks = repmat(struct('nodes',[],'P',[],'skeleton_nodes',[], ...
    'candidate_dofs',[]),nparts,1);
part_patch_counts = accumarray(part_id,1);
part_node_counts = zeros(nparts,1);
local_skeleton_counts = zeros(nparts,1);
local_errors = zeros(nparts,1);
local_norms = zeros(nparts,1);
candidate_nodes = [];
interpolation_entries = 0;
interpolation_bound = 0;

for part = 1:nparts
    patches = find(part_id==part);
    assert(all(conncomp(graph(adjacency(patches,patches)))==1), ...
        'An outgoing skeleton source part is disconnected.')
    nodes = find(node_part==part).';
    source = node_subset(S,nodes);
    w = sqrt(repelem(source.wts(:),3));
    K = proxy_w.*nrccie_matrix(source,T,zk);
    local_norms(part) = norm(K./w.','fro');
    [Jlocal,P,local_errors(part)] = nrccie_outgoing_skeleton(source,K,tol/2);
    selected = nodes((Jlocal(1:3:end)+2)/3);
    blocks(part).nodes = nodes;
    blocks(part).P = P;
    blocks(part).skeleton_nodes = selected;
    blocks(part).candidate_dofs = 3*numel(candidate_nodes)+(1:numel(Jlocal));
    candidate_nodes = [candidate_nodes,selected];
    part_node_counts(part) = numel(nodes);
    local_skeleton_counts(part) = numel(selected);
    interpolation_entries = interpolation_entries+numel(P);
    % In weighted density coordinates this maps a part to its candidates.
    L = (w(Jlocal).*P)./w.';
    interpolation_bound = max(interpolation_bound,norm(L,'fro'));
    clear K L P
end

% Merge candidate columns on the same proxy. Tighten the merge tolerance
% to account for amplification by the block-diagonal local interpolation.
source = node_subset(S,candidate_nodes);
w = sqrt(repelem(source.wts(:),3));
K = proxy_w.*nrccie_matrix(source,T,zk);
candidate_norm = norm(K./w.','fro');
full_norm = norm(local_norms);
merge_tol = min(tol/2,(tol/2)*full_norm/(candidate_norm*interpolation_bound));
[Jmerge,Pmerge,merge_error] = nrccie_outgoing_skeleton(source,K,merge_tol);
skeleton_nodes = candidate_nodes((Jmerge(1:3:end)+2)/3);
Jout = reshape(3*skeleton_nodes+(-2:0).',1,[]);
proxy_error_bound = norm(local_errors.*local_norms)/full_norm ...
    +merge_error*candidate_norm*interpolation_bound/full_norm;
assert(proxy_error_bound<=1.1*tol, ...
    'The composed outgoing proxy error bound exceeds the requested tolerance.')

info.part_id = part_id;
info.adjacency = adjacency;
info.part_patch_counts = part_patch_counts;
info.part_node_counts = part_node_counts;
info.local_skeleton_counts = local_skeleton_counts;
info.candidate_nodes = candidate_nodes;
info.local_errors = local_errors;
info.merge_error = merge_error;
info.merge_tol = merge_tol;
info.proxy_error_bound = proxy_error_bound;
info.largest_proxy_nodes = max([part_node_counts;numel(candidate_nodes)]);
info.full_proxy_mib = 16*(6*size(T.r,2))*(3*S.npts)/2^20;
info.largest_proxy_mib = 16*(6*size(T.r,2))*(3*info.largest_proxy_nodes)/2^20;
info.interpolation_mib = 16*(interpolation_entries+numel(Pmerge))/2^20;
end

function T = node_subset(S,nodes)
T.npts = numel(nodes);
T.r = S.r(:,nodes);
T.dru = S.dru(:,nodes);
T.drv = S.drv(:,nodes);
T.wts = S.wts(nodes);
end
