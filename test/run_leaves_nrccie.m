% NRCCIE self-convergence runs: regenerate the same geometry at each order.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(fullfile(root,'src'))
clear pth dir

%% Parameters
norders = 10;
settings.geometry_options = leaf_geometry_options(); % override shared defaults here
settings.geometry_revision = leaf_geometry_revision();
settings.thickness = .6925; % mm
settings.zk = 2*pi/(10*69.25);
settings.alpha = 1;
settings.eps_quad = 1e-11;
settings.eps_fmm = 1e-9;
settings.eps_gmres = 1e-8;
settings.gmres_restart = [];
settings.gmres_maxit = 1000;
settings.quad_batch_size = 2000;
settings.source_info.r = [0;0;10*69.25];
settings.source_info.edips = -69.25^3*[1;1i;0];

if ~isfolder(fullfile(root,'data')), mkdir(fullfile(root,'data')); end
for norder = norders
    run_order(root,norder,settings);
end

function run_order(root,norder,settings)
%% Same patches and analytic geometry, new quadrature nodes
settings.norder = norder;
opts = settings.geometry_options;
opts.norder = norder;
[S,parts] = generate_leaves_surfer(settings.thickness,opts);
settings.geometry_options = parts.options;
common = cell(1,4); local = cell(1,4);
reflection = cell(1,8);
common_reflection = leaf_midplane_permutation(parts.local_common_surfer,parts.common_groups);
for j = 1:4
    common{j} = affine_transf(parts.local_common_surfer,parts.rotation(:,:,j),parts.shift(:,j));
    local{j} = affine_transf(parts.local_unique_surfers{j},parts.rotation(:,:,j),parts.shift(:,j));
    reflection{2*j-1} = common_reflection;
    reflection{2*j} = leaf_midplane_permutation(parts.local_unique_surfers{j},parts.local_groups{j});
end
settings.rotations = parts.rotation;
settings.shifts = parts.shift;
clear parts
fprintf('Order %d: %d nodes, %d unknowns\n',norder,S.npts,3*S.npts)
zk = settings.zk; alpha = settings.alpha;
eps_quad = settings.eps_quad; eps_fmm = settings.eps_fmm;
eps_gmres = settings.eps_gmres; source_info = settings.source_info;
output_file = fullfile(root,'data',sprintf('run_leaves_nrccie_order%d.mat',norder));

%% Corrections and NRCCIE
timer = tic;
[corrections,correction_info] = common_local_quad_corr_mats( ...
    common,local,eps_quad,zk,settings.quad_batch_size,reflection);
quadrature_time = toc(timer);
clear common local
fprintf('Stored / expanded correction nonzeros: %d / %d\n', ...
    correction_info.stored_nonzeros,correction_info.assembled_nonzeros)
fprintf('Correction storage: %.2f GiB\n',correction_info.stored_bytes/2^30)
operator.npts = S.npts;
operator.r = S.r;
operator.wts = S.wts(:).';
operator.n = S.n;
operator.ru = S.du./vecnorm(S.du,2,1);
operator.zk = zk;
operator.alpha = alpha;
operator.eps_fmm = eps_fmm;
operator.rv = cross(operator.n,operator.ru,1);
operator.apply_corrections = @(d)common_local_apply_quad_corr(corrections,d);
matvec = @(x)apply_nrccie(x,operator);
[einc,hinc] = em3d.incoming_sources(zk,source_info,S,'electric dipole');
normal_einc = sum(S.n.*einc,1);
tangent_rhs = cross(S.n,hinc,1)-alpha*(S.n.*normal_einc-einc);
rhs_components = [sum(operator.ru.*tangent_rhs,1); ...
    sum(operator.rv.*tangent_rhs,1);normal_einc];
rhs = rhs_components(:);

%% Solve and save
fprintf('Order %d: starting GMRES for %d unknowns\n',norder,numel(rhs));
timer = tic;
[solution,flag,relres,iter,resvec] = gmres_with_progress(matvec,rhs,settings.gmres_restart,eps_gmres,settings.gmres_maxit);
solve_time = toc(timer);
fprintf('Order %d: verifying the true residual\n',norder);
true_relative_residual = norm(matvec(solution)-rhs)/norm(rhs);
clear matvec
operator = rmfield(operator,'apply_corrections');
clear corrections
d = reshape(solution,3,S.npts);
surface_current = operator.ru.*d(1,:)+operator.rv.*d(2,:);
surface_charge = d(3,:);
solver_data.gmres_flag = flag;
solver_data.gmres_relative_residual = relres;
solver_data.gmres_iterations = iter;
solver_data.number_of_iterations = numel(resvec)-1;
solver_data.gmres_residual_history = resvec;
solver_data.true_relative_residual = true_relative_residual;
solver_data.quadrature_time = quadrature_time;
solver_data.solve_time = solve_time;
solver_data.correction_info = correction_info;
save(output_file,'S','settings','solver_data','rhs_components', ...
    'surface_current','surface_charge','-v7.3')
fprintf('GMRES flag %d; recomputed residual %.3e; saved %s\n', ...
    flag,true_relative_residual,output_file)
if flag ~= 0 || true_relative_residual > 5*eps_gmres
    error('Order %d GMRES did not meet the requested residual.',norder)
end
end
