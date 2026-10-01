function [J,P,relative_error] = nrccie_incoming_skeleton(S,B,tol)
%NRCCIE_INCOMING_SKELETON Pivoted-QR row ID, completed to whole nodes.
% B approximately equals P*B(J,:).
% B contains incident NRCCIE RHS samples, with three rows per surface node.
% Standard QR pivots select scalar DOFs; their nodes are then completed.
% Each returned node contributes all three consecutive DOFs to J.
% The selected node indices are (J(1:3:end)+2)/3.
% P maps unweighted skeleton values to unweighted values at all DOFs.
% tol controls the area-weighted relative Frobenius error on these samples;
% validate on independent directions before using P for other incident waves.
w = sqrt(repelem(S.wts(:),3));
Bw = w.*B;
sample_norm = norm(Bw,'fro');
if sample_norm==0
    J = [];
    P = zeros(size(B,1),0,'like',B);
    relative_error = 0;
    return
end
[~,R,permutation] = qr(Bw.','econ','vector');
tail_squared = [flipud(cumsum(flipud(sum(abs(R).^2,2))));0];
rank_in = find(tail_squared<=tol^2*tail_squared(1),1)-1;
scalar_dofs = permutation(1:rank_in);
coefficients = R(1:rank_in,1:rank_in)\R(1:rank_in,rank_in+1:end);

% Embed the standard QR interpolant into a skeleton of complete nodes.
nodes = unique(ceil(scalar_dofs/3),'stable');
J = reshape(3*nodes+(-2:0).',1,[]);
[~,columns] = ismember(scalar_dofs,J);
P = zeros(size(B,1),numel(J),'like',B);
P(permutation,columns) = [eye(rank_in);coefficients.'];
P = (P.*w(J).')./w;
% Extra components reproduce themselves; the other rows keep the QR ID.
P(J,:) = eye(numel(J));
relative_error = norm(w.*(B-P*B(J,:)),'fro')/sample_norm;
end
