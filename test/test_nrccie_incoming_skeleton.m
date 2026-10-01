% Small incoming-skeleton test. No integral-equation assembly or solves.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(23)

S = geometries.ellipsoid([1.2,0.85,0.65],[1,1,1],[0;0;0],5,11);
zk = pi/1.2; % One wavelength across the longest diameter.
alpha = 1;
tol = 1e-9;
rotations = sample_so3(512);
validation_rotations = sample_so3(64);
for j = 1:size(rotations,3)
    R = rotations(:,:,j);
    assert(norm(R.'*R-eye(3),'fro')<1e-13 && abs(det(R)-1)<1e-13)
end

% Compare both polarizations against the library's plane-wave evaluator.
Bcheck = nrccie_planewave_rhs(S,zk,alpha,rotations(:,:,1));
d = rotations(:,3,1);
theta = acos(d(3));
phi = atan2(d(2),d(1));
e_theta = [cos(theta)*cos(phi);cos(theta)*sin(phi);-sin(theta)];
e_phi = [-sin(phi);cos(phi);0];
for polarization = 1:2
    p = rotations(:,polarization,1);
    pol = [e_theta.'*p;e_phi.'*p];
    [E,H] = em3d.planewave(zk,[theta;phi],pol,S);
    tangent = cross(S.n,H,1)-alpha*cross(S.n,cross(S.n,E,1),1);
    rhs = [sum(S.dru.*tangent,1);sum(S.drv.*tangent,1);sum(S.n.*E,1)];
    assert(norm(Bcheck(:,polarization)-rhs(:))/norm(rhs(:))<1e-13)
end

% Build a complete-node skeleton and test on fresh incident directions.
Bin = nrccie_planewave_rhs(S,zk,alpha,rotations);
[Jin,Pin,training_error] = nrccie_incoming_skeleton(S,Bin,tol);
Btest = nrccie_planewave_rhs(S,zk,alpha,validation_rotations);
w = sqrt(repelem(S.wts(:),3));
residual = w.*(Btest-Pin*Btest(Jin,:));
validation_error = norm(residual,'fro')/norm(w.*Btest,'fro');
worst_wave_error = max(vecnorm(residual)./vecnorm(w.*Btest));
node_dofs = reshape(Jin,3,[]);
skeleton_nodes = (node_dofs(1,:)+2)/3;

assert(all(mod(node_dofs(1,:),3)==1))
assert(isequal(node_dofs,3*skeleton_nodes+(-2:0).'))
assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
assert(norm(Pin(Jin,:)-eye(numel(Jin)),'fro')<1e-12)
assert(training_error<1.1*tol)
assert(validation_error<3*tol && worst_wave_error<10*tol)
fprintf('Surface: %d nodes, %d DOFs; incoming samples: %d waves\n', ...
    S.npts,3*S.npts,size(Bin,2))
fprintf('Incoming skeleton: %d complete nodes, %d DOFs\n', ...
    numel(skeleton_nodes),numel(Jin))
fprintf('Weighted errors: training %.3e, held-out %.3e, worst wave %.3e\n', ...
    training_error,validation_error,worst_wave_error)

% Scalar QR pivots land at node 1, component 1, then node 2, component 2.
% Completing these nodes must not require inverting their redundant rows.
Ssmall.wts = [1;4;9];
Bsmall = [4,0; .2,.1i; 0,.1; .1,.2; 0,1; .1i,.1; ...
    .15,.1; .1,.15i; .25i,.1];
for nsample = 1:2
    [Jsmall,Psmall,small_error] = nrccie_incoming_skeleton( ...
        Ssmall,Bsmall(:,1:nsample),1e-12);
    assert(isequal(Jsmall,1:3*nsample))
    assert(norm(Psmall(Jsmall,:)-eye(numel(Jsmall)),'fro')<1e-12)
    assert(small_error<1e-12)
end
[Jzero,Pzero,zero_error] = nrccie_incoming_skeleton(Ssmall,zeros(9,2),1e-12);
assert(isempty(Jzero) && isequal(size(Pzero),[9,0]) && zero_error==0)
disp('NRCCIE_INCOMING_SKELETON_PASS')
