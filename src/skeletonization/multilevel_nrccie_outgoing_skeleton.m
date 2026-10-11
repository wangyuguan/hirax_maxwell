function [J,blocks,Pmerge,info] = multilevel_nrccie_outgoing_skeleton(S,T,zk,tol,max_patches)
%MULTILEVEL_NRCCIE_OUTGOING_SKELETON Multilevel ID of outgoing E,H fields.
%
% Syntax:
%   [J,blocks,Pmerge,info] = multilevel_nrccie_outgoing_skeleton( ...
%       S,T,zk,tol,max_patches)
%
% Input arguments:
%   S           - source surfer, including patch connectivity and weights.
%   T           - separated enclosing proxy, with positions r and weights wts.
%   zk          - wave number.
%   tol         - relative source/proxy-weighted Frobenius error tolerance.
%   max_patches - maximum patches per connected finest-level region.
%
% Output arguments:
%   J      - final skeleton DOFs, three consecutive DOFs per selected node.
%   blocks - local nodes, density map P, and candidate_dofs per region.
%   Pmerge - map from finest-level candidate densities to final skeletons.
%   info   - hierarchy, timings and bound on the proxy-sample error.
%
% Region p contributes the following equivalent density on J:
%   Pmerge(:,blocks(p).candidate_dofs) * blocks(p).P * density_p.
% Sum these contributions over regions. Deterministic QR uses the same proxy
% and original source weights at every level; validate at independent targets.

proxy_w = sqrt(repelem(T.wts(:),6));
samplefun = @(source)proxy_w.*nrccie_matrix(source,T,zk);
[J,blocks,Pmerge,info] = multilevel_nrccie_skeleton( ...
    S,samplefun,tol,max_patches,'outgoing');
info.proxy_error_bound = info.relative_error_bound;
info.largest_proxy_nodes = info.largest_sample_nodes;
info.full_proxy_mib = info.full_sample_mib;
info.largest_proxy_mib = info.largest_sample_mib;
end
