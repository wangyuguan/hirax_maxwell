function K = nrccie_matrix(S,T,zk)
%NRCCIE_MATRIX Map [ju;jv;rho] density to E,H at separated targets.
% S is the source surfer. T.r contains target points separated from S.
% K is 6*size(T.r,2)-by-3*S.npts, including source quadrature weights.
% Rows are [Ex;Ey;Ez;Hx;Hy;Hz] at each target. No jump or near correction.
nt = size(T.r,2);
K = complex(zeros(6*nt,3*S.npts));
weights = S.wts(:).';
for i = 1:nt
    displacement = T.r(:,i)-S.r;
    distance = vecnorm(displacement,2,1);
    green = exp(1i*zk*distance).*weights./(4*pi*distance);
    gradient = displacement.*(green.*(1i*zk*distance-1)./distance.^2);
    electric_rows = 6*i-5:6*i-3;
    magnetic_rows = 6*i-2:6*i;
    K(electric_rows,1:3:end) = 1i*zk*green.*S.dru;
    K(electric_rows,2:3:end) = 1i*zk*green.*S.drv;
    K(electric_rows,3:3:end) = -gradient;
    K(magnetic_rows,1:3:end) = cross(gradient,S.dru,1);
    K(magnetic_rows,2:3:end) = cross(gradient,S.drv,1);
end
end
