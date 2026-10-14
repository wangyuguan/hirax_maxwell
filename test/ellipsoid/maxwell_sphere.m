function P = maxwell_sphere(radius,norder)
% MAXWELL_SPHERE Gauss-Legendre x uniform quadrature on a centered sphere.
%
%  Input arguments:
%    * radius: sphere radius
%    * norder: number of polar nodes; use 2*norder azimuthal nodes
%
%  Output arguments:
%    * P: r, n, dru, drv (3 by npts), wts, and npts
%        dru and drv are unit tangents with dru x drv = n.

    j = (1:norder-1).';
    beta = j./sqrt(4*j.^2-1);
    [V,D] = eig(diag(beta,1)+diag(beta,-1));
    [z,ids] = sort(diag(D));
    wz = 2*V(1,ids).^2;
    nphi = 2*norder;
    [phi,z] = meshgrid((0:nphi-1)*2*pi/nphi,z);
    phi = phi(:).';
    z = z(:).';
    s = sqrt(1-z.^2);
    P.n = [s.*cos(phi);s.*sin(phi);z];
    P.r = radius*P.n;
    P.dru = [z.*cos(phi);z.*sin(phi);-s];
    P.drv = [-sin(phi);cos(phi);zeros(size(phi))];
    P.wts = radius^2*repmat(wz,1,nphi)*2*pi/nphi;
    P.npts = numel(z);
end
