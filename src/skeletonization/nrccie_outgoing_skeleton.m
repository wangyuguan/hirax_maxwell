function [J,P,relative_error] = nrccie_outgoing_skeleton(S,K,tol)
%NRCCIE_OUTGOING_SKELETON Pivoted-QR column ID, completed to whole nodes.
% K approximately equals K(:,J)*P. K maps density(:) to outgoing samples;
% its 3*S.npts columns are interleaved [ju;jv;rho], with source weights included.
% For farfield compression, use E,H samples on a closed enclosing proxy
% from nrccie_matrix. Target row weights, if desired, go into K first.
% Standard QR pivots select scalar DOFs; their nodes are then completed.
% The selected node indices are (J(1:3:end)+2)/3.
% P*density(:) gives equivalent skeleton densities, not density(J).
% tol controls ||(K-K(:,J)*P)/W||_F / ||K/W||_F, where W has the repeated
% square-root source area weights. Validate with nrccie_matrix at independent
% exterior targets; proxy error alone does not bound error at other targets.
w = sqrt(repelem(S.wts(:),3));
Kw = K./w.';
sample_norm = norm(Kw,'fro');
if sample_norm==0
    J = [];
    P = zeros(0,size(K,2),'like',K);
    relative_error = 0;
    return
end
[~,R,permutation] = qr(Kw,'econ','vector');
tail_squared = [flipud(cumsum(flipud(sum(abs(R).^2,2))));0];
rank_out = find(tail_squared<=tol^2*tail_squared(1),1)-1;
scalar_dofs = permutation(1:rank_out);
coefficients = R(1:rank_out,1:rank_out)\R(1:rank_out,rank_out+1:end);

% Embed the standard QR interpolant into a skeleton of complete nodes.
nodes = unique(ceil(scalar_dofs/3),'stable');
J = reshape(3*nodes+(-2:0).',1,[]);
[~,rows] = ismember(scalar_dofs,J);
P = zeros(numel(J),size(K,2),'like',K);
P(rows,permutation) = [eye(rank_out),coefficients];
P = (P.*w.')./w(J);
% Extra components reproduce themselves; the other columns keep the QR ID.
P(:,J) = eye(numel(J));
relative_error = norm((K-K(:,J)*P)./w.','fro')/sample_norm;
end
