function R = nrccie_rhs_matrix(S,alpha)
% NRCCIE_RHS_MATRIX Map incident E,H samples to the PEC NRCCIE right hand side.
%
%  Input arguments:
%    * S: physical surfer with n, dru, drv, and npts
%    * alpha: NRCCIE parameter
%
%  Output arguments:
%    * R: sparse 3*npts by 6*npts matrix
%        Input is interleaved [Ex,Ey,Ez,Hx,Hy,Hz]; output is [bu,bv,brho].

    n = S.npts;
    ir = zeros(18,n);
    ic = zeros(18,n);
    values = zeros(18,n);
    for j = 1:n
      u = S.dru(:,j);
      v = S.drv(:,j);
      normal = S.n(:,j);
      block = [alpha*u.',cross(u,normal).';
               alpha*v.',cross(v,normal).';
               normal.',zeros(1,3)];
      [rows,cols] = ndgrid(3*j-2:3*j,6*j-5:6*j);
      ir(:,j) = rows(:);
      ic(:,j) = cols(:);
      values(:,j) = block(:);
    end
    R = sparse(ir(:),ic(:),values(:),3*n,6*n);
end
