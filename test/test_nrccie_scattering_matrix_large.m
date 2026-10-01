% Original 7986-node ellipsoid: manufactured exterior Maxwell scattering map.
% No dense system or global incoming/outgoing interpolation is formed.
% Independent FMM/GMRES separates discretization from compression errors.
% This is a boundary-data-to-radiating-field map, not PEC plane-wave scattering.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir
cache_dir = fullfile(root,'diagnostics','nrccie_scattering_large');
output_dir = fullfile(root,'diagnostics','nrccie_scattering_large_proxy_tol1e8');
if ~isfolder(output_dir), mkdir(output_dir); end
diary(fullfile(output_dir,'run.log'))
diary_cleanup = onCleanup(@()diary('off'));
rng(43)

geometry = struct('semi_axes',[1.2,0.85,0.65],'patches_per_axis',[3,3,4], ...
    'center',[0;0;0],'order',10,'iptype',11);
S = geometries.ellipsoid(geometry.semi_axes,geometry.patches_per_axis, ...
    geometry.center,geometry.order,geometry.iptype);
assert(S.npts==7986)
ndof = 3*S.npts;
zk = 0.9;
alpha = 1;
eps_quad = 1e-11;
factor_tol = 1e-8;
skeleton_tol = 1e-7;
max_nodes = 100;
max_patches = 6;
batch_size = 17;
far_radii = [1e4,2e4];
opts = struct('zk',zk,'alpha',alpha,'jump',-0.5,'fmm',false,'verb',1,'proxy',true);
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = 1e-11;
settings = struct('zk',zk,'alpha',alpha,'jump',opts.jump, ...
    'eps_quad',eps_quad,'eps_fmm',fmm_opts.eps_fmm,'factor_tol',factor_tol, ...
    'skeleton_tol',skeleton_tol,'max_nodes',max_nodes,'max_patches',max_patches, ...
    'batch_size',batch_size,'proxy_radius',4,'proxy_order',6, ...
    'training_positions',128,'training_radius',0.25, ...
    'finite_target_radius',40,'far_radii',far_radii, ...
    'gmres_restart',100,'gmres_tol',1e-10,'gmres_maxit',3,'flam_proxy',true);
fprintf('Original ellipsoid: %d nodes, %d DOFs; dense A would occupy %.2f GiB.\n', ...
    S.npts,ndof,16*ndof^2/1024^3)

% Reuse data independent of the FLAM settings; preserve both no-proxy runs.
cache_settings = rmfield(settings,{'factor_tol','flam_proxy'});
cache_file = fullfile(cache_dir,'quadrature.mat');
reuse_quadrature = false;
if isfile(cache_file)
    cached = load(cache_file,'geometry','settings');
    reuse_quadrature = isequal(cached.geometry,geometry) && ...
        isequal(rmfield(cached.settings,'factor_tol'),cache_settings);
end
if reuse_quadrature
    cached = load(cache_file,'C','quadrature_seconds');
    C = cached.C;
    quadrature_seconds = cached.quadrature_seconds;
    clear cached
    fprintf('Reusing matching quadrature corrections from the 1e-4 run.\n')
else
    fprintf('Building quadrature corrections ...\n')
    quadrature_timer = tic;
    C = nrccie_quad_corr_block(S,eps_quad,zk,[],100);
    quadrature_seconds = toc(quadrature_timer);
end
assert(sum(cellfun(@nnz,C))>0)
save(fullfile(output_dir,'quadrature.mat'),'C','geometry','settings', ...
    'quadrature_seconds','-v7.3')
fprintf('Quadrature corrections: %.2f s; checkpoint saved.\n',quadrature_seconds)

fprintf('Factoring the %d-by-%d exterior NRCCIE system at tolerance %.1e ...\n', ...
    ndof,ndof,factor_tol)
factor_timer = tic;
F = nrccie_rskelf(S,C,opts,max_nodes,factor_tol);
factor_seconds = toc(factor_timer);
assert(F.N==ndof && F.symm=='n')
save(fullfile(output_dir,'factorization.mat'),'F','geometry','settings', ...
    'factor_seconds','-v7.3')
factor_storage = whos('F');
below_root = 0;
for j = 1:F.lvp(end-1)
    below_root = below_root+numel(F.factors(j).rd);
end
fprintf('Factor: %.2f s, %.2f GiB; %d levels, %d DOFs eliminated below root.\n', ...
    factor_seconds,factor_storage.bytes/1024^3,F.nlvl,below_root)
assert(below_root>0,'No low-rank elimination occurred before the root.')

fprintf('Checking rskelf_mv and rskelf_sv against the corrected FMM ...\n')
q = randn(ndof,3)+1i*randn(ndof,3);
q = q./vecnorm(q);
factored_values = rskelf_mv(F,q);
fmm_values = complex(zeros(size(q)));
solution = rskelf_sv(F,q);
solved_values = complex(zeros(size(q)));
for j = 1:size(q,2)
    fmm_values(:,j) = nrccie_apply(S,q(:,j),C,fmm_opts);
    solved_values(:,j) = nrccie_apply(S,solution(:,j),C,fmm_opts);
end
multiply_error = relative_error(factored_values,fmm_values);
solve_residual = relative_error(solved_values,q);
fprintf('rskelf_mv/FMM %.3e; rskelf_sv FMM residual %.3e.\n', ...
    multiply_error,solve_residual)
assert(multiply_error<30*factor_tol && solve_residual<100*factor_tol)
clear q factored_values fmm_values solution solved_values

% Train on outgoing traces of interior dipoles. The callback constructs
% only the requested connected part or candidate union, never the full B.
ntraining = settings.training_positions;
z = 1-2*((0:ntraining-1)+0.5)/ntraining;
phi = pi*(3-sqrt(5))*(0:ntraining-1);
training_positions = settings.training_radius* ...
    [sqrt(1-z.^2).*cos(phi);sqrt(1-z.^2).*sin(phi);z];
rhsfun = @(T)dipole_training_rhs(T,zk,alpha,training_positions);
cache_file = fullfile(cache_dir,'skeletons.mat');
reuse_skeletons = false;
if isfile(cache_file)
    cached = load(cache_file,'geometry','settings','training_positions');
    reuse_skeletons = isequal(cached.geometry,geometry) && ...
        isequal(rmfield(cached.settings,'factor_tol'),cache_settings) && ...
        isequal(cached.training_positions,training_positions);
end
if reuse_skeletons
    load(cache_file,'Jin','Jout','blocks_in','Pmerge_in','blocks_out','Pmerge_out', ...
        'info_in','info_out','incoming_seconds','outgoing_seconds')
    clear cached
    fprintf('Reusing matching skeletons: %d incoming nodes, %d outgoing nodes.\n', ...
        numel(Jin)/3,numel(Jout)/3)
else
    fprintf('Selecting incoming skeleton from %d dipole RHS columns, by parts ...\n', ...
        6*ntraining)
    incoming_timer = tic;
    [Jin,blocks_in,Pmerge_in,info_in] = hierarchical_nrccie_incoming_skeleton( ...
        S,zk,alpha,[],skeleton_tol,max_patches,rhsfun);
    incoming_seconds = toc(incoming_timer);
    fprintf('Incoming: %d nodes from %d connected parts, bound %.3e, %.2f s.\n', ...
        numel(Jin)/3,numel(blocks_in),info_in.training_error_bound,incoming_seconds)

    proxy = geometries.ellipsoid(4*[1,1,1],[1,1,1],[0;0;0],6,11);
    fprintf('Selecting outgoing skeleton with %d proxy nodes ...\n',proxy.npts)
    outgoing_timer = tic;
    [Jout,blocks_out,Pmerge_out,info_out] = hierarchical_nrccie_outgoing_skeleton( ...
        S,proxy,zk,skeleton_tol,max_patches);
    outgoing_seconds = toc(outgoing_timer);
    fprintf('Outgoing: %d nodes from %d connected parts, bound %.3e, %.2f s.\n', ...
        numel(Jout)/3,numel(blocks_out),info_out.proxy_error_bound,outgoing_seconds)
end
assert(isequal(sort([blocks_in.nodes]),1:S.npts))
assert(isequal(sort([blocks_out.nodes]),1:S.npts))
assert(all(info_in.part_patch_counts<=max_patches))
assert(all(info_out.part_patch_counts<=max_patches))
assert(isequal(reshape(Jin,3,[]),3*ceil(Jin(1:3:end)/3)+(-2:0).'))
assert(isequal(reshape(Jout,3,[]),3*ceil(Jout(1:3:end)/3)+(-2:0).'))
save(fullfile(output_dir,'skeletons.mat'),'Jin','Jout','blocks_in','Pmerge_in', ...
    'blocks_out','Pmerge_out','info_in','info_out','training_positions', ...
    'incoming_seconds','outgoing_seconds','geometry','settings','-v7.3')

fprintf('Forming %d-by-%d scattering matrix in batches of %d ...\n', ...
    numel(Jout),numel(Jin),batch_size)
scattering_timer = tic;
Scat = nrccie_scattering_matrix( ...
    F,blocks_in,Pmerge_in,blocks_out,Pmerge_out,batch_size);
scattering_seconds = toc(scattering_timer);
assert(isequal(size(Scat),[numel(Jout),numel(Jin)]))
save(fullfile(output_dir,'scattering_checkpoint.mat'),'Scat','Jin','Jout', ...
    'scattering_seconds','geometry','settings','-v7.3')
fprintf('Scattering matrix formed in %.2f s; checkpoint saved.\n',scattering_seconds)

% A held-out dipole provides an exact radiating exterior solution.
source.r = [0.14;-0.09;0.07];
source.edips = [1;0.25i;-0.3];
source.hdips = [-0.2i;0.35;0.15i];
[boundary_e,boundary_h] = em3d.incoming_sources(zk,source,S,'ehd');
Btest = manufactured_rhs(S,boundary_e,boundary_h,alpha);
Btest = Btest(:);
source_w = sqrt(repelem(S.wts(:),3));
candidate_rhs = Pmerge_in*Btest(Jin);
reconstructed_rhs = complex(zeros(ndof,1));
for part = 1:numel(blocks_in)
    b = blocks_in(part);
    dofs = reshape(3*b.nodes+(-2:0).',1,[]);
    reconstructed_rhs(dofs) = b.P*candidate_rhs(b.candidate_dofs);
end
incoming_error = relative_error(source_w.*reconstructed_rhs,source_w.*Btest);
fprintf('Held-out incoming RHS error: %.3e.\n',incoming_error)
assert(incoming_error<10*skeleton_tol)

fprintf('Solving independent held-out FMM/GMRES reference ...\n')
reference_timer = tic;
Afun = @(x)nrccie_apply(S,x,C,fmm_opts);
[reference_density,gmres_flag,gmres_relres,gmres_iter,gmres_resvec] = ...
    gmres(Afun,Btest,settings.gmres_restart,settings.gmres_tol,settings.gmres_maxit);
reference_seconds = toc(reference_timer);
reference_residual = relative_error(Afun(reference_density),Btest);
factor_density = rskelf_sv(F,Btest);
factor_density_error = relative_error(factor_density,reference_density);
heldout_factor_residual = relative_error(Afun(factor_density),Btest);
fprintf('GMRES flag %d, iterations [%d %d], residual %.3e, %.2f s.\n', ...
    gmres_flag,gmres_iter(1),gmres_iter(2),reference_residual,reference_seconds)
fprintf('Held-out factor density error %.3e; FMM residual %.3e.\n', ...
    factor_density_error,heldout_factor_residual)
save(fullfile(output_dir,'reference.mat'),'reference_density','factor_density', ...
    'Btest','source','gmres_flag','gmres_relres','gmres_iter','gmres_resvec', ...
    'reference_residual','reference_seconds','geometry','settings','-v7.3')
assert(gmres_flag==0 && reference_residual<5*settings.gmres_tol)
assert(heldout_factor_residual<100*factor_tol)

equivalent_density = Scat*Btest(Jin);
reference_equivalent_density = complex(zeros(numel(Jout),1));
for part = 1:numel(blocks_out)
    b = blocks_out(part);
    dofs = reshape(3*b.nodes+(-2:0).',1,[]);
    reference_equivalent_density = reference_equivalent_density+ ...
        Pmerge_out(:,b.candidate_dofs)*(b.P*reference_density(dofs));
end
outgoing_source = node_subset(S,ceil(Jout(1:3:end)/3));
rng(43) % Target directions are reproducible independently of system size.
target_rotations = sample_so3(32);
directions = reshape(target_rotations(:,3,:),3,[]);
evaluation_radii = [settings.finite_target_radius,far_radii];
reference_fields = complex(zeros(6*size(directions,2),numel(evaluation_radii)));
factor_fields = reference_fields;
scattering_fields = reference_fields;
outgoing_fields = reference_fields;
fprintf('Evaluating finite and far fields on 32 independent directions ...\n')
for ir = 1:numel(evaluation_radii)
    targets.r = evaluation_radii(ir)*directions;
    Kout = nrccie_matrix(outgoing_source,targets,zk);
    scattering_fields(:,ir) = Kout*equivalent_density;
    outgoing_fields(:,ir) = Kout*reference_equivalent_density;
    for part = 1:numel(blocks_out)
        nodes = blocks_out(part).nodes;
        dofs = reshape(3*nodes+(-2:0).',1,[]);
        Kpart = nrccie_matrix(node_subset(S,nodes),targets,zk);
        reference_fields(:,ir) = reference_fields(:,ir)+Kpart*reference_density(dofs);
        factor_fields(:,ir) = factor_fields(:,ir)+Kpart*factor_density(dofs);
    end
end
targets.r = settings.finite_target_radius*directions;
[exact_e,exact_h] = em3d.incoming_sources(zk,source,targets,'ehd');
analytic_field = [exact_e;exact_h];
finite_field_error = relative_error(reshape(scattering_fields(:,1),6,[]),analytic_field);
reference_finite_error = relative_error(reshape(reference_fields(:,1),6,[]),analytic_field);
outgoing_field_error = relative_error(outgoing_fields(:,1),reference_fields(:,1));
factor_field_error = relative_error(factor_fields(:,1),reference_fields(:,1));
scattering_field_error = relative_error(scattering_fields(:,1),reference_fields(:,1));

% Remove exp(ikR)/R and extrapolate the leading 1/R correction away.
scale = far_radii.*exp(-1i*zk*far_radii);
far_scattering = scattering_fields(:,2:3).*scale;
far_reference = reference_fields(:,2:3).*scale;
far_factor = factor_fields(:,2:3).*scale;
farfield_scattering = reshape(2*far_scattering(:,2)-far_scattering(:,1),6,[]);
farfield_reference = reshape(2*far_reference(:,2)-far_reference(:,1),6,[]);
farfield_factor = reshape(2*far_factor(:,2)-far_factor(:,1),6,[]);

% Independent dipole asymptotics in the convention of incoming_sources.
phase = exp(-1i*zk*(source.r.'*directions))/(4*pi);
p = repmat(source.edips,1,size(directions,2));
m = repmat(source.hdips,1,size(directions,2));
electric_farfield = -zk^2*phase.*( ...
    p-directions.*sum(directions.*p,1)+cross(directions,m,1));
magnetic_farfield = zk^2*phase.*( ...
    m-directions.*sum(directions.*m,1)-cross(directions,p,1));
farfield_analytic = [electric_farfield;magnetic_farfield];
farfield_e_error = relative_error(farfield_scattering(1:3,:),electric_farfield);
farfield_h_error = relative_error(farfield_scattering(4:6,:),magnetic_farfield);
reference_farfield_e_error = relative_error(farfield_reference(1:3,:),electric_farfield);
reference_farfield_h_error = relative_error(farfield_reference(4:6,:),magnetic_farfield);
farfield_compression_error = relative_error(farfield_scattering,farfield_reference);
farfield_factor_error = relative_error(farfield_factor,farfield_reference);
fprintf('Analytic finite field: scattering %.3e, FMM/GMRES reference %.3e.\n', ...
    finite_field_error,reference_finite_error)
fprintf('Analytic farfield: scattering E %.3e, H %.3e; reference E %.3e, H %.3e.\n', ...
    farfield_e_error,farfield_h_error,reference_farfield_e_error,reference_farfield_h_error)
fprintf('Finite-field errors vs reference: outgoing %.3e, factor %.3e, Scat %.3e.\n', ...
    outgoing_field_error,factor_field_error,scattering_field_error)
fprintf('Farfield errors vs reference: factor %.3e, Scat %.3e.\n', ...
    farfield_factor_error,farfield_compression_error)

errors = struct('factor_multiply',multiply_error,'factor_residual',solve_residual, ...
    'incoming_training_bound',info_in.training_error_bound, ...
    'outgoing_proxy_bound',info_out.proxy_error_bound,'incoming_heldout',incoming_error, ...
    'reference_residual',reference_residual,'factor_density',factor_density_error, ...
    'heldout_factor_residual',heldout_factor_residual,'outgoing_field',outgoing_field_error, ...
    'factor_field',factor_field_error,'scattering_field',scattering_field_error, ...
    'finite_field',finite_field_error,'reference_finite_field',reference_finite_error, ...
    'farfield_e',farfield_e_error,'farfield_h',farfield_h_error, ...
    'reference_farfield_e',reference_farfield_e_error, ...
    'reference_farfield_h',reference_farfield_h_error, ...
    'farfield_factor',farfield_factor_error,'farfield_compression',farfield_compression_error);
stats = struct('surface_nodes',S.npts,'factor_levels',F.nlvl, ...
    'eliminated_below_root',below_root,'factor_bytes',factor_storage.bytes, ...
    'dense_bytes',16*ndof^2,'quadrature_seconds',quadrature_seconds, ...
    'factor_seconds',factor_seconds,'incoming_seconds',incoming_seconds, ...
    'outgoing_seconds',outgoing_seconds,'scattering_seconds',scattering_seconds, ...
    'reference_seconds',reference_seconds,'incoming_nodes',numel(Jin)/3, ...
    'outgoing_nodes',numel(Jout)/3);
save(fullfile(output_dir,'scattering_results.mat'), ...
    'Scat','Jin','Jout','geometry','settings','errors','stats', ...
    'blocks_in','Pmerge_in','blocks_out','Pmerge_out','info_in','info_out', ...
    'source','training_positions','directions','farfield_analytic', ...
    'farfield_reference','farfield_factor','farfield_scattering', ...
    'analytic_field','reference_fields','scattering_fields', ...
    'gmres_flag','gmres_relres','gmres_iter','gmres_resvec','-v7.3')
assert(outgoing_field_error<200*skeleton_tol)
assert(scattering_field_error<100*factor_tol+200*skeleton_tol)
assert(farfield_compression_error<100*factor_tol+200*skeleton_tol)
assert(max([reference_finite_error,reference_farfield_e_error, ...
    reference_farfield_h_error])<1e-5,'Fine-grid analytic reference is inaccurate.')
assert(max([finite_field_error,farfield_e_error,farfield_h_error])<100*factor_tol)
disp('NRCCIE_SCATTERING_MATRIX_LARGE_PASS')
diary off


function B = dipole_training_rhs(T,zk,alpha,positions)
B = complex(zeros(3*T.npts,6*size(positions,2)));
moments = eye(3);
for j = 1:size(positions,2)
    source.r = positions(:,j);
    for component = 1:6
        if component<=3
            source.edips = moments(:,component);
            type = 'ed';
        else
            source.hdips = moments(:,component-3);
            type = 'hd';
        end
        [E,H] = em3d.incoming_sources(zk,source,T,type);
        rhs = manufactured_rhs(T,E,H,alpha);
        B(:,6*(j-1)+component) = rhs(:);
    end
end
end


function rhs = manufactured_rhs(S,E,H,alpha)
normal_e = sum(S.n.*E,1);
tangent = -cross(S.n,H,1)+alpha*(S.n.*normal_e-E);
rhs = [sum(S.dru.*tangent,1);sum(S.drv.*tangent,1);-normal_e];
end


function T = node_subset(S,nodes)
T.npts = numel(nodes);
T.r = S.r(:,nodes);
T.dru = S.dru(:,nodes);
T.drv = S.drv(:,nodes);
T.wts = S.wts(nodes);
end


function err = relative_error(value,reference)
err = norm(value-reference,'fro')/max(norm(reference,'fro'),realmin);
end
