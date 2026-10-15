clear all
clc

run(fullfile(fileparts(mfilename('fullpath')),'startup.m'))


narms = 5;
amp = .5;
wavelength = 1.0;
zk = 2*pi/wavelength;


unscaled_starfish_diameter = ...
    2*(1+amp)*cos(pi/(2*narms));
starfish_scale = wavelength/unscaled_starfish_diameter;

starfish_cparams = struct( ...
    'nchmin',40,'maxchunklen',4/zk,'eps',1e-12);
starfish_pref = struct('k',16);
starfish_chnkr = chunkerfunc( ...
    @(t) starfish(t,narms,amp,[0;0],0,starfish_scale), ...
    starfish_cparams,starfish_pref);

proxy_radius = 1.25*wavelength;
proxy_cparams = struct( ...
    'nchmin',16,'maxchunklen',4/zk,'eps',1e-12);
circle_chnkr = chunkerfunc( ...
    @(t) ellipse(t,proxy_radius,proxy_radius),proxy_cparams);

fprintf('starfish discretization: %d chunks, %d points\n', ...
    starfish_chnkr.nch,starfish_chnkr.npt);
fprintf('wavelength = %.3f, wavenumber = %.3f\n',wavelength,zk);
fprintf('approximate starfish diameter = %.3f wavelengths\n', ...
    unscaled_starfish_diameter*starfish_scale/wavelength);


figure(1)
plot(starfish_chnkr)
hold on
plot(circle_chnkr)
axis equal
legend('starfish','proxy circle','Location','best')


% Sound-soft Helmholtz scattering with the robust CFIE D - i*k*S.
cfie_coefficients = [1,-1i*zk];
helm_cfie = kernel('helmholtz','c',zk,cfie_coefficients);


opts = struct('forcesmooth',true);
Kout = chunkerkernevalmat( ...
    starfish_chnkr,helm_cfie,circle_chnkr,opts);


% Solve the exterior Dirichlet BIE for an incident plane wave

nstar = starfish_chnkr.npt;
xstar = starfish_chnkr.r(:,:);
incident_direction = [1;0.35];
incident_direction = incident_direction/norm(incident_direction);
wavevector = zk*incident_direction;
uinc_boundary = exp(1i*(wavevector.'*xstar)).';

CFIE_boundary = chunkermat(starfish_chnkr,helm_cfie);
A = 0.5*eye(nstar) + CFIE_boundary;
mu = A\(-uinc_boundary);

boundary_residual = norm( ...
    A*mu + uinc_boundary) ...
    /norm(uinc_boundary);

fprintf('BIE boundary residual       = %.3e\n',boundary_residual);


%% Find incoming skeleton points from uniformly sampled plane waves

id_tolerance = 1e-10;
nplane = 128;
plane_angles = 2*pi*(0:nplane-1)/nplane;
plane_directions = [cos(plane_angles);sin(plane_angles)];

% Bin(j,m) is the trace at boundary point j of the plane wave arriving
% from direction m.  A column ID of Bin' is a row ID of Bin, so Jin
% contains physical points on the starfish boundary.
Bin = exp(1i*zk*(xstar.'*plane_directions));
[Jin,Jin_redundant,Tin] = id(Bin',id_tolerance);

rank_in = numel(Jin);
Pin = zeros(nstar,rank_in,'like',Bin);
Pin(Jin,:) = eye(rank_in);
Pin(Jin_redundant,:) = Tin';

Bin_from_skeleton = Pin*Bin(Jin,:);
incoming_training_error = norm( ...
    Bin-Bin_from_skeleton,'fro')/norm(Bin,'fro');

% Test directions halfway between the directions used to construct the ID.
test_plane_angles = plane_angles + pi/nplane;
test_plane_directions = ...
    [cos(test_plane_angles);sin(test_plane_angles)];
Bin_test = exp(1i*zk*(xstar.'*test_plane_directions));
Bin_test_from_skeleton = Pin*Bin_test(Jin,:);
incoming_test_error = norm( ...
    Bin_test-Bin_test_from_skeleton,'fro')/norm(Bin_test,'fro');

% The incident wave used by the BIE above was not put on the angle grid.
uinc_from_skeleton = Pin*uinc_boundary(Jin);
incoming_wave_error = norm( ...
    uinc_boundary-uinc_from_skeleton)/norm(uinc_boundary);

fprintf('incoming skeleton rank      = %d of %d\n',rank_in,nstar);
fprintf('incoming training error     = %.3e\n',incoming_training_error);
fprintf('incoming off-grid error     = %.3e\n',incoming_test_error);
fprintf('current incident-wave error = %.3e\n',incoming_wave_error);


%% Find outgoing skeleton points by a column ID of Kout

[~,R,permutation] = qr(Kout,'vector');

relative_R_diagonal = abs(diag(R))/abs(R(1,1));
rank_id = sum(relative_R_diagonal > id_tolerance);
Jout = permutation(1:rank_id);

R11 = R(1:rank_id,1:rank_id);
R12 = R(1:rank_id,rank_id+1:end);
T = R11\R12;

Pout = zeros(rank_id,size(Kout,2));
Pout(:,permutation) = [eye(rank_id),T];


%% Compare the full and skeletonized outgoing fields

equivalent_strength = Pout*mu;

u_circle_full = Kout*mu;
u_circle_skeleton = ...
    Kout(:,Jout)*equivalent_strength;

relative_field_error = norm( ...
    u_circle_full-u_circle_skeleton) ...
    /norm(u_circle_full);

relative_operator_error = norm( ...
    Kout-Kout(:,Jout)*Pout,'fro') ...
    /norm(Kout,'fro');

fprintf('outgoing skeleton rank      = %d of %d\n',rank_id,nstar);
fprintf('operator ID relative error  = %.3e\n',relative_operator_error);
fprintf('field relative error        = %.3e\n',relative_field_error);


%% Plot the skeleton points and the field comparison

figure(2)
clf
plot(starfish_chnkr)
hold on
scatter(xstar(1,Jout),xstar(2,Jout),60,'filled')
axis equal
legend('starfish','outgoing skeleton points','Location','best')
title(sprintf('Outgoing skeleton points: %d of %d',rank_id,nstar))

figure(4)
clf
plot(starfish_chnkr)
hold on
scatter(xstar(1,Jin),xstar(2,Jin),60,'filled')
axis equal
legend('starfish','incoming skeleton points','Location','best')
title(sprintf('Incoming skeleton points: %d of %d',rank_in,nstar))

xcircle = circle_chnkr.r(:,:);
theta_circle = mod(atan2(xcircle(2,:),xcircle(1,:)),2*pi);
[theta_circle,angle_order] = sort(theta_circle);

figure(3)
clf
subplot(2,1,1)
plot(theta_circle,real(u_circle_full(angle_order)),'LineWidth',1.5)
hold on
plot(theta_circle,real(u_circle_skeleton(angle_order)), ...
    '--','LineWidth',1.5)
xlabel('angle on proxy circle')
ylabel('real part of scattered field')
legend('full','skeleton approximation','Location','best')
grid on

subplot(2,1,2)
semilogy(theta_circle, ...
    abs(u_circle_full(angle_order)-u_circle_skeleton(angle_order)), ...
    'LineWidth',1.5)
xlabel('angle on proxy circle')
ylabel('absolute error')
grid on


[Q,R,p] = qr(Kout,'vector');

relative_diagonal = abs(diag(R))/abs(R(1,1));
tol = 1e-8;


figure
semilogy(relative_diagonal,'o-')
yline(tol,'--')
xlabel('Pivot number')
ylabel('|R_{jj}|/|R_{11}|')
grid on



r = sum(relative_diagonal > tol);
Jout = p(1:r);


xstar = starfish_chnkr.r(:,:);
figure
plot(starfish_chnkr)
hold on
scatter(xstar(1,Jout),xstar(2,Jout), 60,'filled')
axis equal




% R11 = R(1:r,1:r);
% R12 = R(1:r,r+1:end);
%
% T = R11\R12;
%
% Pout = zeros(r,size(Kout,2));
% Pout(:,p) = [eye(r),T];
