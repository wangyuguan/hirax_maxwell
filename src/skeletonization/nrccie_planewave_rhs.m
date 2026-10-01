function B = nrccie_planewave_rhs(S,zk,alpha,rotations)
%NRCCIE_PLANEWAVE_RHS PEC scattering RHS for sampled incident plane waves.
% R(:,3) is the propagation direction; R(:,1:2) are the E polarizations.
% E = p exp(ik d.x), H = (d x p) exp(ik d.x), with normalized impedance.
% B is 3*S.npts-by-2*size(rotations,3), interleaved by boundary node.
% Columns 2*j-1 and 2*j are the two polarizations of rotations(:,:,j).
% Uses the physical PEC RHS paired with nrccie_apply's jump=+0.5.
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
