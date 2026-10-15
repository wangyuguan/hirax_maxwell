function [J,P,relative_error] = nrccie_incoming_skeleton(S,B,tol)
%NRCCIE_INCOMING_SKELETON Pivoted-QR row ID, completed to whole nodes.
%
% Input arguments:
%   S - source surfer with positive quadrature weights.
%   B - incident NRCCIE RHS samples, with three rows per surface node.
%   tol - relative area-weighted Frobenius error tolerance.
%
% Output arguments:
%   J - selected DOFs, including all three components of each selected node.
%   P - interpolation matrix, B approximately equals P*B(J,:).
%   relative_error - achieved weighted error on the supplied samples.
%
% P maps unweighted values. Validate on independent incident waves.
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
