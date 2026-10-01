function F = nrccie_rskelf(S,C,opts,max_nodes,tol)
%NRCCIE_RSKELF Factor the corrected self NRCCIE system with FLAM.
% C is {Cslp,Cx,Cy,Cz}; opts supplies zk, alpha, and jump; optional verb=1
% prints FLAM's per-level progress. A local proxy is used by default;
% opts.proxy=false restores sampling of all actual off-diagonal entries.
% opts.proxy_order optionally sets the proxy's polar sampling order.
% max_nodes is the leaf occupancy target in surface nodes (3 DOFs each).
% tol is FLAM's local ID tolerance, not a guaranteed global solve tolerance.
% The proxy covers both directions of the nonsymmetric operator. Actual
% matrix blocks still include quadrature corrections and the identity jump.
% FLAM's internal pivots are scalar DOFs; external incoming/outgoing skeletons
% still contain complete nodes.
% Use rskelf_mv(F,x) for A*x and rskelf_sv(F,b) for A\b.
matrix_block = @(rows,cols)nrccie_system_block(S,[],C,opts,rows,cols);
coordinates = repelem(S.r,1,3);
flam_options = struct('symm','n');
if isfield(opts,'verb')
    flam_options.verb = opts.verb;
end
proxyfun = [];
if ~isfield(opts,'proxy') || opts.proxy
    proxyfun = nrccie_rskelf_proxy(S,C,opts,tol);
end
F = rskelf(matrix_block,coordinates,3*max_nodes,tol,proxyfun,flam_options);
end
