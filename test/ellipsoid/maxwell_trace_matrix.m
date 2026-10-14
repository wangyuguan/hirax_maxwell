function L = maxwell_trace_matrix(P)
% MAXWELL_TRACE_MATRIX Extract tangential electric and magnetic fields.
%
%  Input arguments:
%    * P: surface with unit tangent frames dru, drv and npts
%
%  Output arguments:
%    * L: sparse 4*npts by 6*npts matrix
%        L maps interleaved [Ex,Ey,Ez,Hx,Hy,Hz] to [Eu,Ev,Hu,Hv].

    n = P.npts;
    nodes = 0:n-1;
    rows = reshape((1:4*n).',4,n);
    ir = repmat(rows,3,1);
    ic = [6*nodes+1;6*nodes+1;6*nodes+4;6*nodes+4;
          6*nodes+2;6*nodes+2;6*nodes+5;6*nodes+5;
          6*nodes+3;6*nodes+3;6*nodes+6;6*nodes+6];
    values = [P.dru(1,:);P.drv(1,:);P.dru(1,:);P.drv(1,:);
              P.dru(2,:);P.drv(2,:);P.dru(2,:);P.drv(2,:);
              P.dru(3,:);P.drv(3,:);P.dru(3,:);P.drv(3,:)];
    L = sparse(ir(:),ic(:),values(:),4*n,6*n);
end
