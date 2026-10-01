function [Jin,blocks,Pmerge,info] = hierarchical_nrccie_incoming_skeleton(S,zk,alpha,rotations,tol,max_patches,rhsfun)
%HIERARCHICAL_NRCCIE_INCOMING_SKELETON Connected parts with shared RHS samples.
% Each part and its candidate union use standard QR row ID, completed to
% whole nodes. The same RHS samples are used throughout, and original surface
% weights are preserved. Jin contains global interleaved RHS DOFs.
% For RHS values b_final at the final skeleton, reconstruct each part with
%   bcandidate = Pmerge*b_final;
%   for i = 1:numel(blocks)
%       b_i = blocks(i).P*bcandidate(blocks(i).candidate_dofs,:);
%   end
% RHS matrices are built for one part or the candidate union at a time;
% the global interpolation is kept as local factors and Pmerge.
% info.training_error_bound bounds the area-weighted Frobenius error on the
% sampled data; independently validate accuracy on new boundary data.
% Optional rhsfun(T) supplies a 3*T.npts-by-nsample matrix on a node subset.
% Its columns must represent the same data throughout. For this mode pass
% rotations=[]; otherwise the default samples the supplied plane waves.

if nargin<7 || isempty(rhsfun)
    rhsfun = @(T)nrccie_planewave_rhs(T,zk,alpha,rotations);
end

[part_id,adjacency] = partition_surface_patches(S,max_patches);
nparts = max(part_id);
node_part = part_id(S.patch_id);
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
        'An incoming skeleton surface part is disconnected.')
    nodes = find(node_part==part).';
    source = node_subset(S,nodes);
    w = sqrt(repelem(source.wts(:),3));
    B = rhsfun(source);
    local_norms(part) = norm(w.*B,'fro');
    [Jlocal,P,local_errors(part)] = nrccie_incoming_skeleton(source,B,tol/2);
    selected = nodes((Jlocal(1:3:end)+2)/3);
    blocks(part).nodes = nodes;
    blocks(part).P = P;
    blocks(part).skeleton_nodes = selected;
    blocks(part).candidate_dofs = 3*numel(candidate_nodes)+(1:numel(Jlocal));
    candidate_nodes = [candidate_nodes,selected];
    part_node_counts(part) = numel(nodes);
    local_skeleton_counts(part) = numel(selected);
    interpolation_entries = interpolation_entries+numel(P);
    % In weighted RHS coordinates this maps candidates to their full part.
    L = (w.*P)./w(Jlocal).';
    interpolation_bound = max(interpolation_bound,norm(L,'fro'));
    clear B L P
end

% Merge rows at candidate nodes. Tighten the merge tolerance to account for
% amplification by the block-diagonal local interpolation.
source = node_subset(S,candidate_nodes);
w = sqrt(repelem(source.wts(:),3));
B = rhsfun(source);
candidate_norm = norm(w.*B,'fro');
full_norm = norm(local_norms);
merge_tol = min(tol/2,(tol/2)*full_norm/(candidate_norm*interpolation_bound));
[Jmerge,Pmerge,merge_error] = nrccie_incoming_skeleton(source,B,merge_tol);
skeleton_nodes = candidate_nodes((Jmerge(1:3:end)+2)/3);
Jin = reshape(3*skeleton_nodes+(-2:0).',1,[]);
training_error_bound = norm(local_errors.*local_norms)/full_norm ...
    +merge_error*candidate_norm*interpolation_bound/full_norm;
assert(training_error_bound<=1.1*tol, ...
    'The composed incoming training error bound exceeds the requested tolerance.')

info.part_id = part_id;
info.adjacency = adjacency;
info.part_patch_counts = part_patch_counts;
info.part_node_counts = part_node_counts;
info.local_skeleton_counts = local_skeleton_counts;
info.candidate_nodes = candidate_nodes;
info.local_errors = local_errors;
info.merge_error = merge_error;
info.merge_tol = merge_tol;
info.training_error_bound = training_error_bound;
info.largest_rhs_nodes = max([part_node_counts;numel(candidate_nodes)]);
info.full_rhs_mib = 16*(3*S.npts)*size(B,2)/2^20;
info.largest_rhs_mib = 16*(3*info.largest_rhs_nodes)*size(B,2)/2^20;
info.interpolation_mib = 16*(interpolation_entries+numel(Pmerge))/2^20;
end

function T = node_subset(S,nodes)
T.npts = numel(nodes);
T.r = S.r(:,nodes);
T.dru = S.dru(:,nodes);
T.drv = S.drv(:,nodes);
T.n = S.n(:,nodes);
T.wts = S.wts(nodes);
end
