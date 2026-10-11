% Four-leaf NRCCIE scattering: direct RS-S of the complete corrected A.
% Edit the parameters below, then run this script in MATLAB.
clear
test_dir = fileparts(mfilename('fullpath'));
run(fullfile(test_dir,'..','startup.m'))

%% Parameters
tol = 1e-8;                       % one tolerance throughout
norder = 4;
leaf_radius = 69.25;              % mm (original HIRAX scale)
thickness = 0.6925;               % mm
wavelength = 138.5;               % one nominal leaf width, mm
zk = 2*pi/wavelength;
alpha = 1;
jump = 0.5;
run_reference = false;            % optional plain full-system GMRES
geometry_only = false;

max_nodes = 256;                  % RS-S leaf-box occupancy in physical nodes
tree_axes = [1,2];                % thin geometry; operator remains 3D
max_refinement = 3;               % direct residual corrections; no GMRES
skeleton_max_patches = 16;
ntraining = 1024;                 % directions, two polarizations each
nheldout = 2;                     % four held-out polarized RHSs
nfield = 48;
proxy_order = 10;                 % 726 proxy points
proxy_radius = 4*wavelength;
field_radii = [7,8,10]*wavelength;
quad_batch_size = 128;
scattering_batch_size = 24;
maxiter = 2000;                   % total GMRES iterations; no restarts

output_dir = fullfile(test_dir,'results',sprintf('four_leaf_order%d_%s', ...
    norder,char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))));
for name = {'factor.mat','result.mat','geometry.mat'}
    assert(~isfile(fullfile(output_dir,name{1})),'Use a fresh output directory.')
end
if ~isfolder(output_dir), mkdir(output_dir); end
run_timer = tic;
timings = struct(); errors = struct(); stats = struct();
factor_file = fullfile(output_dir,'factor.mat');
result_file = fullfile(output_dir,'result.mat');
fprintf('Solver: direct RS-S of complete A; source: %s\n',test_dir)

%% Geometry
% One master leaf, rotated into four copies with matched nodes and frames.
timer = tic;
geometry_opts.norder = norder;
geometry_opts.rim_width = 0.028128271246*leaf_radius;
geometry_opts.cap_collar_width = geometry_opts.rim_width;
geometry_opts.cap_mesh_spacing = 0.5*leaf_radius;
geometry_opts.cap_mesh_spacing_center = geometry_opts.cap_mesh_spacing;
geometry_opts.cap_mesh_spacing_side = geometry_opts.cap_mesh_spacing;
geometry_opts.cap_mesh_side_start = 0.62;
geometry_opts.outline_refinement = 1;
geometry_opts.wall_profile_refinement = 1;
geometry_opts.surface_gap = 4;
geometry_opts.chunkie_order = 20;
geometry_opts.chunkie_n0 = 3;
geometry_opts.chunkie_nchs = 3;
geometry_opts.chunkie_newton_iterations = 30;
[S,parts] = generate_leaves_surfer(thickness,geometry_opts);
leaves = parts.leaves;
geometry = struct('kind','leaf','order',norder,'cycle',[1,2,4,3], ...
    'base',parts.base_surfer,'rotations',parts.rotation,'shifts',parts.shift, ...
    'require_cross_near',true,'geometry_revision',parts.geometry_revision, ...
    'options',parts.options,'thickness',thickness,'geometry_seconds',toc(timer));
clear parts
timings.geometry = geometry.geometry_seconds;
if geometry_only
    save(fullfile(output_dir,'geometry.mat'),'S','leaves','geometry','-v7.3')
    disp('NRCCIE_SCATTERING_MATRIX_LEAF_GEOMETRY_PASS'); return
end
timer = tic;
errors.geometry = check_geometry(S,leaves,geometry.cycle);
center = (min(S.r,[],2)+max(S.r,[],2))/2;
extent = max(vecnorm(S.r-center,2,1));
assert(proxy_radius>extent && all(field_radii>proxy_radius))
N = 3*S.npts;
stats.nodes_per_leaf = leaves{1}.npts; stats.total_dofs = N;
stats.reference_enabled = run_reference;
timings.geometry_checks = toc(timer);
fprintf('Four %s surfaces: %d nodes/leaf, %d DOFs, tolerance %.1e.\n', ...
    geometry.kind,stats.nodes_per_leaf,N,tol)

%% Quadrature and full-system operator
% Reuse C4 near corrections; the reference operator has no low-rank factors.
timer = tic;
[C4,quad_info] = fourfold_quad_corr_mats(leaves{1},leaves{2},leaves{3},leaves{4}, ...
    tol,zk,quad_batch_size);
timings.c4_quadrature = toc(timer);
timings.c4_quadrature_blocks = quad_info.block_seconds;
assert(all(sum(quad_info.base_block_nonzeros([2,4],:),2)>0))
storage = whos('C4'); stats.quadrature_bytes = storage.bytes;
% C4 quadrature reuse is always used for both the solver and reference.
opts.zk = zk;
opts.alpha = alpha;
opts.jump = jump;
opts.fmm = false;
opts.verb = 1;
opts.tree_axes = tree_axes;
fmm_opts = opts; fmm_opts.fmm = true; fmm_opts.eps_fmm = tol/100;
Cref = @(rho)fourfold_apply_quad_corr(C4,rho);
Afun = @(rho)nrccie_apply(S,rho,Cref,fmm_opts);
stats.reference_quadrature_bytes = 0;
errors.quadrature = []; errors.sampled_blocks = [];

%% Factorization: direct RS-S of the complete corrected four-leaf A
checkpoint('factor_started',output_dir,stats,errors,timings);
timer = tic;
[F,factor_info] = nrccie_srskelf(S,leaves,C4,opts,max_nodes,tol);
timings.factor = toc(timer);
[stats,timings] = save_factor( ...
    F,factor_info,opts,tol,geometry,stats,timings,factor_file);
clear C4 quad_info
index = (1:N).'; rhs = sin(sqrt(2)*index)+1i*cos(sqrt(3)*index); rhs = rhs/norm(rhs);
[~,solve_info] = nrccie_direct_solve(F,rhs,Afun,tol,max_refinement);
timings.test_rhs_solve = solve_info.seconds;
timings.test_rhs_backsolve = solve_info.backsolve_seconds;
timings.test_rhs_refinement = solve_info.refinement_seconds;
timings.test_rhs_residual = solve_info.residual_seconds;
stats.test_rhs_refinement = solve_info.refinement_steps;
errors.test_rhs_initial_residual = solve_info.initial_true_relres;
errors.test_rhs_residual = solve_info.true_relres;
fprintf('Test RHS: %d refinements, initial/final residual %.3e/%.3e.\n', ...
    solve_info.refinement_steps,solve_info.initial_true_relres, ...
    solve_info.true_relres)
assert(errors.test_rhs_residual<=max(tol,100*eps), ...
    'Test-RHS true residual is too large.')
clear index rhs

%% Incoming skeleton: one leaf, reused four times
checkpoint('incoming_started',output_dir,stats,errors,timings);
timer = tic;
training_rotations = nrccie_deterministic_rotations(ntraining,0,true);
rhsfun = @(T)nrccie_planewave_rhs(T,zk,alpha,training_rotations);
[Jin,blocks_in,Pmerge_in,info_in] = fourfold_nrccie_skeleton( ...
    leaves,rhsfun,tol,skeleton_max_patches,'incoming');
timings.incoming = toc(timer);

%% Outgoing skeleton: fields on the proxy sphere
timer = tic;
proxy = geometries.sphere(proxy_radius,1,center,proxy_order,11);
assert(min(vecnorm(proxy.r-center,2,1))>extent)
errors.proxy_covariance = check_proxy_covariance(proxy,leaves,zk);
proxy_w = sqrt(repelem(proxy.wts(:),6));
samplefun = @(T)proxy_w.*nrccie_matrix(T,proxy,zk);
[Jout,blocks_out,Pmerge_out,info_out] = fourfold_nrccie_skeleton( ...
    leaves,samplefun,tol,skeleton_max_patches,'outgoing');
timings.outgoing = toc(timer);
Sin = node_subset(S,ceil(Jin(1:3:end)/3));
Sout = node_subset(S,ceil(Jout(1:3:end)/3));
[stats,errors] = skeleton_metadata(stats,errors,info_in,info_out);

%% Fresh plane waves and the full scattering matrix
timer = tic;
heldout_rotations = nrccie_deterministic_rotations(nheldout,0.371);
Btest = nrccie_planewave_rhs(S,zk,alpha,heldout_rotations);
candidate_rhs = Pmerge_in*Btest(Jin,:); reconstructed = complex(zeros(size(Btest)));
for p = 1:numel(blocks_in)
    block = blocks_in(p); dofs = reshape(3*block.nodes+(-2:0).',1,[]);
    reconstructed(dofs,:) = block.P*candidate_rhs(block.candidate_dofs,:);
end
w = sqrt(repelem(S.wts(:),3));
errors.incoming_heldout = column_errors(w.*reconstructed,w.*Btest);
clear reconstructed candidate_rhs
timings.incoming_heldout = toc(timer);
fprintf('Incoming held-out maximum %.3e.\n',max(errors.incoming_heldout))
assert(max(errors.incoming_heldout)<10*tol,'Incoming held-out reconstruction failed.')
checkpoint('scattering_started',output_dir,stats,errors,timings);
timer = tic;
% Scat = Pout * A^{-1} * Pin, using direct RS-S solves in batches.
[Scat,scat_info] = nrccie_scattering_matrix( ...
    F,Afun,tol,max_refinement,blocks_in,Pmerge_in,blocks_out,Pmerge_out, ...
    scattering_batch_size);
timings.scattering = toc(timer);
assert(isequal(size(Scat),[numel(Jout),numel(Jin)]))
stats.scattering_size = size(Scat);
stats.scattering_solver = 'direct-a-rss';
stats.scattering_max_refinement = max(cellfun(@(b)max(b.refinement_steps),scat_info.batches));
timings.scattering_backsolve = sum(cellfun(@(b)b.backsolve_seconds,scat_info.batches));
timings.scattering_refinement = sum(cellfun(@(b)b.refinement_seconds,scat_info.batches));
timings.scattering_residual = sum(cellfun(@(b)b.residual_seconds,scat_info.batches));
errors.scattering_initial_residual = max(cellfun( ...
    @(b)max(b.initial_true_relres),scat_info.batches));
errors.scattering_residual = max(cellfun(@(b)max(b.true_relres),scat_info.batches));
clear scat_info

%% Solve four held-out RHSs; optionally compare with plain full-system GMRES
checkpoint('heldout_solves_started',output_dir,stats,errors,timings);
nwave = size(Btest,2); assert(nwave==4)
reference_density = [];
stats.gmres_flag = []; stats.gmres_iter = []; errors.reference_residual = [];
timings.reference_solve_per_wave = []; timings.reference_residual_per_wave = [];
[direct_density,direct_info] = nrccie_direct_solve( ...
    F,Btest,Afun,tol,max_refinement);
timings.direct_solve = direct_info.seconds;
timings.direct_backsolve = direct_info.backsolve_seconds;
timings.direct_refinement = direct_info.refinement_seconds;
timings.direct_residual = direct_info.residual_seconds;
stats.direct_refinement = direct_info.refinement_steps;
errors.direct_initial_residual = direct_info.initial_true_relres;
errors.direct_residual = direct_info.true_relres;
fprintf(['Direct A RS-S: refinement %d/%d, initial/final max residual ' ...
    '%.3e/%.3e, %.2f s.\n'],min(direct_info.refinement_steps), ...
    max(direct_info.refinement_steps),max(direct_info.initial_true_relres), ...
    max(direct_info.true_relres),direct_info.seconds)
if run_reference
    reference_density = complex(zeros(N,nwave));
    stats.gmres_flag = zeros(1,nwave); stats.gmres_iter = zeros(nwave,2);
end
for wave = 1:nwave
    if run_reference
        timer = tic;
        [reference_density(:,wave),stats.gmres_flag(wave),reported,stats.gmres_iter(wave,:)] = ...
            gmres(Afun,Btest(:,wave),[],tol,min(maxiter,N),[],[],zeros(N,1));
        timings.reference_solve_per_wave(wave) = toc(timer);
        fprintf('Plain GMRES RHS %d: flag %d, residual %.3e, %.2f s.\n', ...
            wave,stats.gmres_flag(wave),reported,timings.reference_solve_per_wave(wave))
        timer = tic;
        errors.reference_residual(wave) = relative_error(Afun(reference_density(:,wave)),Btest(:,wave));
        timings.reference_residual_per_wave(wave) = toc(timer);
    end
end
timings.reference_solve = []; errors.density = []; timings.density_comparison = [];
if run_reference
    timings.reference_solve = sum(timings.reference_solve_per_wave);
    timer = tic; errors.density = column_errors(direct_density,reference_density);
    timings.density_comparison = toc(timer);
end
timings.residual_checks = sum(timings.reference_residual_per_wave)+ ...
    sum(timings.density_comparison);

%% Actual E/H fields at 7, 8 and 10 wavelengths, away from training targets
timer = tic;
field_rotations = nrccie_deterministic_rotations(nfield,0.613);
directions = reshape(field_rotations(:,3,:),3,[]);
equivalent_density = Scat*Btest(Jin,:);
[errors.fields,fields] = nrccie_scattering_field_checks(S,Sout,zk, ...
    center,directions,field_radii,reference_density,direct_density,equivalent_density);
timings.fields = toc(timer);
timings.field_radii_seconds = errors.fields.radius_seconds;

%% Save first, then mark successful only after all numerical checks pass
run_info = struct('matlab',version,'finished',char(datetime('now')), ...
    'reference_enabled',run_reference,'external_skeletons','one leaf, C4 reuse', ...
    'scope','same-discretization solver/compression test; not mesh convergence');
result = struct('passed',false,'kind',geometry.kind,'errors',errors,'stats',stats, ...
    'timings',timings,'factor_file',factor_file,'result_file',result_file);
surface_nodes = S.r;
timer = tic;
save(result_file,'Scat','Sin','Sout','Jin','Jout','surface_nodes','tol','zk','alpha','opts','run_info', ...
    'norder','thickness','geometry_opts','max_nodes','skeleton_max_patches', ...
    'ntraining','nheldout','nfield','proxy_order','quad_batch_size', ...
    'scattering_batch_size','max_refinement','maxiter', ...
    'center','wavelength','proxy_radius','field_radii','heldout_rotations','fields','result','-v7.3')
timings.result_save = toc(timer); timings.total = toc(run_timer);
result.timings = timings; save(result_file,'result','-append')
solve_limit = max(tol,100*eps); field_limit = max(100*tol,100*eps);
reference_limit = max(5*tol,50*eps);
assert(errors.test_rhs_residual<=solve_limit,'Test-RHS true residual is too large.')
assert(max(errors.incoming_heldout)<10*tol,'Incoming held-out reconstruction failed.')
assert(max(errors.direct_residual)<=solve_limit,'Direct true residual is too large.')
assert(errors.scattering_residual<=solve_limit,'Scat solves did not meet the A residual target.')
f = errors.fields;
assert(max([f.direct_vs_scattering_e(:);f.direct_vs_scattering_h(:)])<field_limit, ...
    'Scat and direct finite-radius E/H fields disagree.')
if run_reference
    assert(all(stats.gmres_flag==0) && max(errors.reference_residual)<reference_limit, ...
        'Plain full-system GMRES did not converge accurately.')
    assert(max(errors.density)<field_limit,'Direct and reference densities disagree.')
    assert(max([f.reference_vs_direct_e(:);f.reference_vs_direct_h(:); ...
        f.reference_vs_scattering_e(:);f.reference_vs_scattering_h(:)])<field_limit, ...
        'Finite-radius E/H fields disagree with plain full-system GMRES.')
end
result.passed = true; save(result_file,'result','-append')
checkpoint('passed',output_dir,stats,errors,timings);
print_timings(timings);
fprintf('Saved %s\nNRCCIE_SCATTERING_MATRIX_%s_PASS\n',result_file,upper(geometry.kind))


function error = check_geometry(S,leaves,cycle)
assert(iscell(leaves) && numel(leaves)==4 && isequal(cycle,[1,2,4,3]))
assert(S.npts==4*leaves{1}.npts && all(cellfun(@(L)L.npts,leaves)==leaves{1}.npts))
Q = [0,1,0;-1,0,0;0,0,1]; error = 0;
for p = 1:4
    a = leaves{cycle(p)}; b = leaves{cycle(mod(p,4)+1)};
    error = max(error,relative_error([b.r,b.n,b.dru,b.drv],Q*[a.r,a.n,a.dru,a.drv]));
end
assert(error<1e-12,'The four surfaces must be rotated copies.')
end


function [stats,timings] = save_factor( ...
        F,factor_info,opts,tol,geometry,stats,timings,filename)
assert(F.N==stats.total_dofs && ~factor_info.randomized && ...
    ~factor_info.proxy && ~factor_info.full_matrix_materialized)
stats.solver = 'direct-a-rss';
stats.factor_bytes = factor_info.factor_bytes;
stats.factor_levels = factor_info.levels;
stats.factor_root_dofs = factor_info.root_dofs;
stats.factor_remaining_by_level = factor_info.remaining_by_level;
stats.factor_local_tolerance = factor_info.local_tolerance;
factor_metadata = struct('kind',geometry.kind,'order',geometry.order, ...
    'geometry_options',geometry.options,'thickness',geometry.thickness, ...
    'geometry_revision',geometry.geometry_revision,'total_dofs',F.N, ...
    'nodes_per_leaf',stats.nodes_per_leaf,'cycle',geometry.cycle, ...
    'factor',factor_info,'seconds',timings.factor);
timer = tic; save(filename,'F','factor_info','opts','tol','factor_metadata','-v7.3')
timings.factor_save = toc(timer);
fprintf('Factor: %.2f s, %.3f GiB, root %d/%d DOFs.\n', ...
    timings.factor,stats.factor_bytes/2^30,stats.factor_root_dofs,F.N)
end


function [stats,errors] = skeleton_metadata(stats,errors,in,out)
assert(~in.randomized && ~out.randomized && in.reused && out.reused)
errors.incoming_training = in.relative_error_bound; errors.outgoing_proxy = out.relative_error_bound;
stats.incoming_nodes_per_leaf = in.skeleton_nodes_per_leaf;
stats.outgoing_nodes_per_leaf = out.skeleton_nodes_per_leaf;
stats.incoming_nodes = sum(in.skeleton_nodes_per_leaf);
stats.outgoing_nodes = sum(out.skeleton_nodes_per_leaf); stats.skeleton_reuse = 'c4';
stats.incoming_sample_mib = in.largest_sample_mib; stats.outgoing_sample_mib = out.largest_sample_mib;
stats.incoming_levels = in.levels; stats.outgoing_levels = out.levels;
fprintf('Per leaf: %d incoming / %d outgoing out of %d nodes; bounds %.3e / %.3e.\n', ...
    stats.incoming_nodes_per_leaf(1),stats.outgoing_nodes_per_leaf(1),stats.nodes_per_leaf, ...
    errors.incoming_training,errors.outgoing_proxy)
end


function error = check_proxy_covariance(proxy,leaves,zk)
Q = [0,1,0;-1,0,0;0,0,1]; rotated = Q*proxy.r;
distance = zeros(1,proxy.npts); permutation = zeros(1,proxy.npts);
for j = 1:proxy.npts
    [distance(j),permutation(j)] = min(vecnorm(proxy.r-rotated(:,j),2,1));
end
scale = max(1,max(abs(proxy.r),[],'all'));
assert(max(distance)<1e3*eps(scale) && numel(unique(permutation))==proxy.npts)
assert(relative_error(proxy.wts(permutation),proxy.wts)<1e-12)
nodes = unique(round(linspace(1,leaves{1}.npts,7)));
K1 = nrccie_matrix(node_subset(leaves{1},nodes),proxy,zk);
K2 = nrccie_matrix(node_subset(leaves{2},nodes),struct('r',proxy.r(:,permutation)),zk);
error = relative_error(K2,kron(speye(proxy.npts),blkdiag(Q,Q))*K1);
assert(error<1e-12,'Outgoing samples must respect C4 covariance.')
end


function source = node_subset(S,nodes)
source = struct('npts',numel(nodes),'r',S.r(:,nodes),'n',S.n(:,nodes), ...
    'dru',S.dru(:,nodes),'drv',S.drv(:,nodes),'wts',S.wts(nodes));
end


function values = column_errors(A,B)
values = vecnorm(A-B,2,1)./max(vecnorm(B,2,1),realmin);
end


function value = relative_error(A,B)
value = norm(A-B,'fro')/max(norm(B,'fro'),realmin);
end


function checkpoint(stage,output_dir,stats,errors,timings)
progress = struct('stage',stage,'passed',strcmp(stage,'passed'), ...
    'updated',char(datetime('now')), ...
    'stats',stats,'errors',errors,'timings',timings);
save(fullfile(output_dir,'progress.mat'),'progress')
fprintf('Stage %s.\n',stage)
end


function print_timings(timings)
fprintf('\nTiming (seconds; details are included in their parent stages)\n')
names = fieldnames(timings);
for j = 1:numel(names)
    name = names{j}; value = timings.(name);
    if isnumeric(value), fprintf('  %-30s %s\n',name,mat2str(value,6)); end
end
end
