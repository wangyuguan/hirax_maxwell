function Q = stratton_chu_matrix(P,T,zk)
% STRATTON_CHU_MATRIX Evaluate fields from tangential traces on a surface.
%
%  Input arguments:
%    * P: closed source surface with r, dru, drv, wts, and npts
%        Unit tangents must satisfy dru x drv = outward normal.
%    * T: separated targets, with positions T.r
%    * zk: nonzero wave number, e^{-i omega t}, epsilon=mu=1
%
%  Output arguments:
%    * Q: 6*nt by 4*P.npts matrix including source weights
%        Input is interleaved [Eu,Ev,Hu,Hv]; output is [Ex,Ey,Ez,Hx,Hy,Hz].
%
%  With J=n x H and M=-n x E, the representation is
%    E = ik (S J + grad div S J / k^2) - curl S M,
%    H = curl S J + ik (S M + grad div S M / k^2).
%  Q reproduces outgoing fields outside P. For a regular incident field,
%  -Q reproduces that field inside P. No on-surface jump is implemented.

    nt = size(T.r,2);
    Q = complex(zeros(6*nt,4*P.npts));
    for i = 1:nt
      displacement = T.r(:,i)-P.r;
      r = vecnorm(displacement,2,1);
      direction = displacement./r;
      green = exp(1i*zk*r).*P.wts./(4*pi*r);
      gradient = direction.*(green.*(1i*zk-1./r));
      a = 1i*zk*green.*(1+1i./(zk*r)-1./(zk*r).^2);
      b = 1i*zk*green.*(-1-3i./(zk*r)+3./(zk*r).^2);
      du = a.*P.dru+b.*direction.*sum(direction.*P.dru,1);
      dv = a.*P.drv+b.*direction.*sum(direction.*P.drv,1);
      cu = cross(gradient,P.dru,1);
      cv = cross(gradient,P.drv,1);
      erows = 6*i-5:6*i-3;
      hrows = 6*i-2:6*i;
      % Eu -> M=-v, Ev -> M=u, Hu -> J=v, Hv -> J=-u.
      Q(erows,1:4:end) = cv;
      Q(erows,2:4:end) = -cu;
      Q(erows,3:4:end) = dv;
      Q(erows,4:4:end) = -du;
      Q(hrows,1:4:end) = -dv;
      Q(hrows,2:4:end) = du;
      Q(hrows,3:4:end) = cv;
      Q(hrows,4:4:end) = -cu;
    end
end
