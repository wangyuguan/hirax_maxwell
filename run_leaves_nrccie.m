% NRCCIE self-convergence runs: regenerate the same geometry at each order.
clear
root=fileparts(mfilename('fullpath'));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'),fullfile(root,'src'))
clear pth dir

%% Parameters
norders=[4 6 8 10 12];
settings.geometry_options=struct(); % override geometry defaults here
settings.thickness=.6925; % mm
settings.zk=2*pi/(10*69.25);
settings.alpha=1;
settings.eps_quad=1e-11;
settings.eps_fmm=1e-9;
settings.eps_gmres=1e-8;
settings.gmres_restart=[]; % unrestarted GMRES
settings.gmres_maxit=1000; % maximum iterations
settings.quad_batch_size=2000;
settings.source_info.r=[0;0;10*69.25];
settings.source_info.edips=-69.25^3*[1;1i;0];

if ~isfolder(fullfile(root,'data')), mkdir(fullfile(root,'data')); end
for norder=norders
    run_order(root,norder,settings);
end

function run_order(root,norder,settings)
%% Same patches and analytic geometry, new quadrature nodes
settings.norder=norder;
opts=settings.geometry_options;
opts.norder=norder;
[S,parts]=hirax_common_local_leaves_surfer(settings.thickness,opts);
settings.geometry_options=parts.options;
common=cell(1,4); local=cell(1,4);
reflection=cell(1,8);
common_reflection=leaf_midplane_permutation(parts.local_common_surfer,parts.common_groups);
for j=1:4
    common{j}=affine_transf(parts.local_common_surfer,parts.rotation(:,:,j),parts.shift(:,j));
    local{j}=affine_transf(parts.local_unique_surfers{j},parts.rotation(:,:,j),parts.shift(:,j));
    reflection{2*j-1}=common_reflection;
    reflection{2*j}=leaf_midplane_permutation(parts.local_unique_surfers{j},parts.local_groups{j});
end
settings.rotations=parts.rotation;
settings.shifts=parts.shift;
clear parts
fprintf('Order %d: %d nodes, %d unknowns\n',norder,S.npts,3*S.npts)
zk=settings.zk; alpha=settings.alpha;
eps_quad=settings.eps_quad; eps_fmm=settings.eps_fmm;
eps_gmres=settings.eps_gmres; source_info=settings.source_info;
output_file=fullfile(root,'data',sprintf('run_leaves_nrccie_order%d.mat',norder));

%% Corrections and NRCCIE
timer=tic;
[corrections,correction_info]=common_local_quad_corr_mats( ...
    common,local,eps_quad,zk,settings.quad_batch_size,reflection);
quadrature_time=toc(timer);
clear common local
fprintf('Stored / expanded correction nonzeros: %d / %d\n', ...
    correction_info.stored_nonzeros,correction_info.assembled_nonzeros)
fprintf('Correction storage: %.2f GiB\n',correction_info.stored_bytes/2^30)
operator=struct('npts',S.npts,'r',S.r,'wts',S.wts(:).', ...
    'n',S.n,'ru',S.du./vecnorm(S.du,2,1), ...
    'zk',zk,'alpha',alpha,'eps_fmm',eps_fmm);
operator.rv=cross(operator.n,operator.ru,1);
operator.apply_corrections=@(d)common_local_apply_quad_corr(corrections,d);
matvec=@(x)apply_nrccie(x,operator);
[einc,hinc]=em3d.incoming_sources(zk,source_info,S,'electric dipole');
normal_einc=sum(S.n.*einc,1);
tangent_rhs=cross(S.n,hinc,1)-alpha*(S.n.*normal_einc-einc);
rhs_components=[sum(operator.ru.*tangent_rhs,1); ...
    sum(operator.rv.*tangent_rhs,1);normal_einc];
rhs=rhs_components(:);

%% Solve and save
timer=tic;
[solution,flag,relres,iter,resvec]=gmres(matvec,rhs,settings.gmres_restart,eps_gmres,settings.gmres_maxit);
solve_time=toc(timer);
true_relative_residual=norm(matvec(solution)-rhs)/norm(rhs);
clear matvec
operator=rmfield(operator,'apply_corrections');
clear corrections
d=reshape(solution,3,S.npts);
surface_current=operator.ru.*d(1,:)+operator.rv.*d(2,:);
surface_charge=d(3,:);
solver_data=struct('gmres_flag',flag,'gmres_relative_residual',relres, ...
    'gmres_iterations',iter,'number_of_iterations',numel(resvec)-1, ...
    'gmres_residual_history',resvec, ...
    'true_relative_residual',true_relative_residual, ...
    'quadrature_time',quadrature_time,'solve_time',solve_time, ...
    'correction_info',correction_info);
save(output_file,'S','settings','solver_data','rhs_components', ...
    'surface_current','surface_charge','-v7.3')
fprintf('GMRES flag %d; recomputed residual %.3e; saved %s\n', ...
    flag,true_relative_residual,output_file)
assert(flag==0 && true_relative_residual<=5*eps_gmres, ...
    'Solve did not converge; inspect this order before continuing.')
end
