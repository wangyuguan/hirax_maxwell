% Incoming skeleton on the ellipsoid from test_nrccie_ellipsoid_analytic.m.
% Only sample and compress incident RHS values; no integral-equation solve.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(23)

% Same physical ellipsoid, wavenumber, and alpha as the analytic test.
semi_axes = [1.2,0.85,0.65];
zk = 0.9;
alpha = 1;
tol = 1e-9;
nsample = 512;             % Two plane-wave polarizations per rotation.
nvalidation = 64;
use_analytic_mesh = false; % Keep local tests small; true uses 7,986 nodes.

if use_analytic_mesh
    patches_per_axis = [3,3,4];
    surface_order = 10;
    mesh_name = 'original analytic-test mesh';
else
    patches_per_axis = [1,1,1];
    surface_order = 5;
    mesh_name = 'small local mesh (same physical ellipsoid)';
end
S = geometries.ellipsoid( ...
    semi_axes,patches_per_axis,[0;0;0],surface_order,11);
fprintf('Mesh: %s\n',mesh_name)
fprintf('Semi-axes [%.2f %.2f %.2f], k %.2f, alpha %.2f, tolerance %.1e\n', ...
    semi_axes,zk,alpha,tol)
fprintf('Surface: %d patches, %d nodes, %d DOFs\n', ...
    S.npatches,S.npts,3*S.npts)

rotations = sample_so3(nsample);
validation_rotations = sample_so3(nvalidation);
Bin = nrccie_planewave_rhs(S,zk,alpha,rotations);
skeleton_timer = tic;
[Jin,Pin,training_error] = nrccie_incoming_skeleton(S,Bin,tol);
skeleton_seconds = toc(skeleton_timer);
skeleton_nodes = (Jin(1:3:end)+2)/3;
skeleton_r = S.r(:,skeleton_nodes);
assert(isequal(reshape(Jin,3,[]),3*skeleton_nodes+(-2:0).'))
assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
assert(norm(Pin(Jin,:)-eye(numel(Jin)),'fro')<1e-12)
assert(training_error<1.1*tol)
clear Bin

% Check interpolation on independent incident directions and polarizations.
Btest = nrccie_planewave_rhs(S,zk,alpha,validation_rotations);
w = sqrt(repelem(S.wts(:),3));
residual = w.*(Btest-Pin*Btest(Jin,:));
validation_error = norm(residual,'fro')/norm(w.*Btest,'fro');
worst_wave_error = max(vecnorm(residual)./vecnorm(w.*Btest));
assert(validation_error<3*tol && worst_wave_error<10*tol)

fprintf('Incoming samples: %d training waves, %d independent test waves\n', ...
    2*nsample,2*nvalidation)
fprintf('Retained nodes: %d / %d (%.2f%%), node compression %.2fx\n', ...
    numel(skeleton_nodes),S.npts,100*numel(skeleton_nodes)/S.npts, ...
    S.npts/numel(skeleton_nodes))
fprintf('Retained DOFs: %d / %d; Pin size: %d x %d\n', ...
    numel(Jin),3*S.npts,size(Pin,1),size(Pin,2))
fprintf('Weighted errors: training %.3e, held-out %.3e, worst wave %.3e\n', ...
    training_error,validation_error,worst_wave_error)
fprintf('Skeleton construction: %.2f s\n',skeleton_seconds)
disp('NRCCIE_ELLIPSOID_SKELETON_PASS')
