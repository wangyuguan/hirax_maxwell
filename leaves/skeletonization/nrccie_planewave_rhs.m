function B = nrccie_planewave_rhs(S,zk,alpha,rotations)
%NRCCIE_PLANEWAVE_RHS PEC scattering RHS for sampled incident plane waves.
%
% Input arguments:
%   S - source surfer.
%   zk, alpha - wave number and NRCCIE coupling parameter.
%   rotations - incident frames, with propagation direction in column 3
%               and electric polarizations in columns 1:2.
%
% Output arguments:
%   B - 3*S.npts-by-2*size(rotations,3) RHS matrix, interleaved by node.
%       Columns 2*j-1 and 2*j use the two polarizations of rotations(:,:,j).
%
% E = p exp(ik d.x), H = (d x p) exp(ik d.x), with normalized impedance.
% The PEC RHS is paired with nrccie_apply's jump=+0.5.
    nsample = size(rotations,3);
    B = complex(zeros(3*S.npts,2*nsample));
    for j = 1:nsample
      R = rotations(:,:,j);
      phase = exp(1i*zk*(R(:,3).'*S.r));
      for polarization = 1:2
        p = R(:,polarization);
        E = p.*phase;
        H = cross(R(:,3),p).*phase;
        normal_e = sum(S.n.*E,1);
        tangent = cross(S.n,H,1)-alpha*(S.n.*normal_e-E);
        rhs = [sum(S.dru.*tangent,1); ...
          sum(S.drv.*tangent,1);normal_e];
        B(:,2*j-2+polarization) = rhs(:);
      end
    end
end
