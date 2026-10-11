function [J,blocks,Pmerge,info] = multilevel_nrccie_incoming_skeleton(S,zk,alpha,rotations,tol,max_patches,rhsfun)
%MULTILEVEL_NRCCIE_INCOMING_SKELETON Multilevel ID of incident boundary data.
%
% Syntax:
%   [J,blocks,Pmerge,info] = multilevel_nrccie_incoming_skeleton( ...
%       S,zk,alpha,rotations,tol,max_patches,rhsfun)
%
% Input arguments:
%   S           - source surfer, including patch connectivity and weights.
%   zk, alpha   - NRCCIE wave number and coupling parameter.
%   rotations   - 3-by-3-by-nwave incident frames, two polarizations each.
%   tol         - relative area-weighted Frobenius sample-error tolerance.
%   max_patches - maximum patches per connected finest-level region.
%   rhsfun      - optional @(T) returning 3*T.npts-by-nsample boundary data;
%                 defaults to the plane waves above. Use the same samples
%                 on every node subset; rotations may be [] with rhsfun.
%
% Output arguments:
%   J      - final skeleton DOFs, three consecutive DOFs per selected node.
%   blocks - local nodes, interpolation P, and candidate_dofs per region.
%   Pmerge - map from final skeleton values to finest-level candidate values.
%   info   - hierarchy, timings and bound on the supplied sample error.
%
% Reconstruct region p from b(J,:) using
%   blocks(p).P * Pmerge(blocks(p).candidate_dofs,:) * b(J,:).
% Every merge uses deterministic QR on its children's skeleton nodes.
% Independently validate new incident data not used in construction.

if nargin<7 || isempty(rhsfun)
    rhsfun = @(T)nrccie_planewave_rhs(T,zk,alpha,rotations);
end
[J,blocks,Pmerge,info] = multilevel_nrccie_skeleton( ...
    S,rhsfun,tol,max_patches,'incoming');
info.training_error_bound = info.relative_error_bound;
info.largest_rhs_nodes = info.largest_sample_nodes;
info.full_rhs_mib = info.full_sample_mib;
info.largest_rhs_mib = info.largest_sample_mib;
end
