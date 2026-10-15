function P = maxwell_ellipsoid(axes,norder,center)
% MAXWELL_ELLIPSOID Gauss-Legendre x uniform quadrature on an ellipsoid.
%
%  Input arguments:
%    * axes: three positive semiaxes, aligned with x, y, and z
%    * norder: number of polar nodes; use 2*norder azimuthal nodes
%    * center: center of the ellipsoid (default [0;0;0])
%
%  Output arguments:
%    * P: r, n, dru, drv (3 by npts), wts, and npts
%        dru and drv are orthonormal unit tangents with dru x drv = n.
%        Traces use the same Eu,Ev,Hu,Hv ordering as MAXWELL_SPHERE.
%
%  The affine image x=center+D*s of the unit sphere has area element
%    dA = det(D)*norm(D^{-T}*s)*dA_s.
%  Surface weights include this factor. The two transformed sphere tangents
%  are generally not orthogonal, so we normalize the first and obtain the
%  second from n x dru instead of using both affine tangents directly.

    if nargin<3 || isempty(center), center = zeros(3,1); end
    validateattributes(axes,{'numeric'},{'real','finite','positive','numel',3});
    validateattributes(norder,{'numeric'},{'scalar','integer','>=',2});
    validateattributes(center,{'numeric'},{'real','finite','numel',3});
    axes = axes(:);
    center = center(:);
    P = maxwell_sphere(1,norder);
    normal = P.n./axes;
    normal_length = vecnorm(normal,2,1);
    first_tangent = axes.*P.dru;
    P.r = center+axes.*P.r;
    P.n = normal./normal_length;
    P.dru = first_tangent./vecnorm(first_tangent,2,1);
    P.drv = cross(P.n,P.dru,1);
    P.wts = prod(axes)*normal_length.*P.wts;
end
