% Small ellipsoid scattering map with a manufactured exterior Maxwell field.
% Report discretization error against analytic dipoles separately from the
% compression error against the same discrete system solved by dense LU.
% This manufactured boundary-data map is not a PEC plane-wave scattering test.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(43)

% Legendre nodes do not coincide across neighboring patches.
S = geometries.ellipsoid([1.2,0.85,0.65],[2,2,2],[0;0;0],4,11);
assert(S.npts==600)
ndof = 3*S.npts;
zk = 0.9;
alpha = 1;
eps_quad = 1e-9;
factor_tol = 1e-4;
skeleton_tol = 1e-7;
max_nodes = 100;
batch_size = 17;
opts = struct('zk',zk,'alpha',alpha,'jump',-0.5,'fmm',false);
fmm_opts = opts;
fmm_opts.fmm = true;
fmm_opts.eps_fmm = 1e-11;
fprintf('Building quadrature corrections for %d ellipsoid nodes ...\n',S.npts)
C = nrccie_quad_corr_block(S,eps_quad,zk,[],100);
assert(sum(cellfun(@nnz,C))>0)
[~,A] = nrccie_apply(S,zeros(ndof,1),C,opts);

% Arbitrary scalar indices: repeated, reordered, and mixed components.
rows = [1,2,3,38,4,ndof-1,2,183];
cols = [2,1,3,39,6,ndof,1,182];
Ablock = nrccie_system_block(S,[],C,opts,rows,cols);
block_error = relative_error(Ablock,A(rows,cols));
assert(block_error<1e-13)
positive_opts = opts;
positive_opts.jump = 0.5;
Aplus = nrccie_system_block(S,[],C,positive_opts,rows,cols);
assert(norm(Aplus-(Ablock+(rows(:)==cols)),'fro')<1e-13)
% Same-node off-diagonal components must not receive an identity jump.
assert(norm(nrccie_system_block(S,[],C,positive_opts,1:3,1:3) ...
    -nrccie_system_block(S,[],C,opts,1:3,1:3)-eye(3),'fro')<1e-13)

% A few nearby exterior targets also exercise rectangular corrections.
target_nodes = [1,8,14,24,39,77,159];
T.r = 1.1*S.r(:,target_nodes);
T.n = S.n(:,target_nodes);
T.dru = S.dru(:,target_nodes);
T.drv = S.drv(:,target_nodes);
T.npts = numel(target_nodes);
Ccross = nrccie_quad_corr_block(S,eps_quad,zk,T);
[~,Across] = nrccie_apply(S,zeros(ndof,1),T,Ccross,opts);
cross_rows = [3*T.npts,1,5,2,5];
cross_block = nrccie_system_block(S,T,Ccross,opts,cross_rows,cols);
assert(relative_error(cross_block,Across(cross_rows,cols))<1e-13)
assert(norm(cross_block-nrccie_system_block( ...
    S,T,Ccross,positive_opts,cross_rows,cols),'fro')==0)

fprintf('Factoring the %d-by-%d exterior NRCCIE system ...\n',ndof,ndof)
factor_timer = tic;
F = nrccie_rskelf(S,C,opts,max_nodes,factor_tol);
factor_seconds = toc(factor_timer);
assert(F.N==ndof && F.symm=='n')
q = randn(ndof,3)+1i*randn(ndof,3);
q = q./vecnorm(q);
factored_values = rskelf_mv(F,q);
fmm_values = complex(zeros(size(q)));
for j = 1:size(q,2)
    fmm_values(:,j) = nrccie_apply(S,q(:,j),C,fmm_opts);
end
fmm_dense_error = relative_error(fmm_values,A*q);
multiply_error = relative_error(factored_values,fmm_values);
assert(fmm_dense_error<5e-8)
assert(multiply_error<30*factor_tol)

solution = rskelf_sv(F,q);
dense_solution = A\q;
solved_values = complex(zeros(size(q)));
for j = 1:size(q,2)
    solved_values(:,j) = nrccie_apply(S,solution(:,j),C,fmm_opts);
end
solve_error = relative_error(solution,dense_solution);
solve_residual = relative_error(solved_values,q);
assert(solve_error<100*factor_tol && solve_residual<100*factor_tol)
cross_fmm = nrccie_apply(S,q(:,1),T,Ccross,fmm_opts);
assert(relative_error(cross_fmm,Across*q(:,1))<5e-8)

% The manufactured data are outgoing traces of interior dipoles. Train on
% these traces, rather than assuming regular plane-wave RHS span this space.
ntraining = 128;
z = 1-2*((0:ntraining-1)+0.5)/ntraining;
phi = pi*(3-sqrt(5))*(0:ntraining-1);
training_positions = 0.25*[sqrt(1-z.^2).*cos(phi);sqrt(1-z.^2).*sin(phi);z];
fprintf('Training incoming skeleton on %d interior-dipole RHS columns ...\n', ...
    6*size(training_positions,2))
B = complex(zeros(ndof,6*size(training_positions,2)));
moments = eye(3);
for j = 1:size(training_positions,2)
    training_source.r = training_positions(:,j);
    for component = 1:6
        if component<=3
            training_source.edips = moments(:,component);
            type = 'ed';
        else
            training_source.hdips = moments(:,component-3);
            type = 'hd';
        end
        [E,H] = em3d.incoming_sources(zk,training_source,S,type);
        rhs = manufactured_rhs(S,E,H,alpha);
        B(:,6*(j-1)+component) = rhs(:);
    end
end
[Jin,Pin,incoming_training_error] = nrccie_incoming_skeleton(S,B,skeleton_tol);
blocks_in = struct('nodes',1:S.npts,'P',Pin, ...
    'skeleton_nodes',ceil(Jin(1:3:end)/3),'candidate_dofs',1:numel(Jin));
Pmerge_in = eye(numel(Jin));
assert(incoming_training_error<1.1*skeleton_tol)

% Outgoing proxy compression still acts on arbitrary surface densities.
% Singleton patches are connected and avoid using coarse extrapolated
% Legendre corners to infer adjacency in this small regression.
proxy = geometries.ellipsoid(4*[1,1,1],[1,1,1],[0;0;0],6,11);
fprintf('Selecting outgoing skeleton with %d proxy nodes ...\n',proxy.npts)
[Jout,blocks_out,Pmerge_out,info_out] = hierarchical_nrccie_outgoing_skeleton( ...
    S,proxy,zk,skeleton_tol,1);
assert(isequal(sort([blocks_in.nodes]),1:S.npts))
assert(isequal(sort([blocks_out.nodes]),1:S.npts))
assert(all(info_out.part_patch_counts==1))
assert(isequal(reshape(Jin,3,[]),3*ceil(Jin(1:3:end)/3)+(-2:0).'))
assert(isequal(reshape(Jout,3,[]),3*ceil(Jout(1:3:end)/3)+(-2:0).'))

% Assemble the full outgoing reference map only for this small test.
Pout = complex(zeros(numel(Jout),ndof));
for part = 1:numel(blocks_out)
    b = blocks_out(part);
    dofs = reshape(3*b.nodes+(-2:0).',1,[]);
    Pout(:,dofs) = Pmerge_out(:,b.candidate_dofs)*b.P;
end
assert(norm(Pin(Jin,:)-eye(numel(Jin)),'fro')<1e-11)
assert(norm(Pout(:,Jout)-eye(numel(Jout)),'fro')<1e-11)
% Always exercise an incomplete final RHS batch, regardless of selected rank.
while mod(numel(Jin),batch_size)==0, batch_size = batch_size+1; end
fprintf('Forming the %d-by-%d scattering matrix in batches of %d ...\n', ...
    numel(Jout),numel(Jin),batch_size)
scattering_timer = tic;
Scat = nrccie_scattering_matrix( ...
    F,blocks_in,Pmerge_in,blocks_out,Pmerge_out,batch_size);
scattering_seconds = toc(scattering_timer);
Scat_reference = Pout*(A\Pin);
scattering_error = relative_error(Scat,Scat_reference);
assert(isequal(size(Scat),[numel(Jout),numel(Jin)]))
assert(scattering_error<100*factor_tol)

% Held-out interior dipole generates an exact radiating exterior solution.
source.r = [0.14;-0.09;0.07];
source.edips = [1;0.25i;-0.3];
source.hdips = [-0.2i;0.35;0.15i];
[boundary_e,boundary_h] = em3d.incoming_sources(zk,source,S,'ehd');
Btest = manufactured_rhs(S,boundary_e,boundary_h,alpha);
Btest = Btest(:);
source_w = sqrt(repelem(S.wts(:),3));
incoming_error = relative_error(source_w.*(Pin*Btest(Jin,:)),source_w.*Btest);
dense_density = A\Btest;
equivalent_density = Scat*Btest(Jin,:);
target_rotations = sample_so3(32);
directions = reshape(target_rotations(:,3,:),3,[]);
targets.r = 40*directions;
Ktest = nrccie_matrix(S,targets,zk);
exact_discrete_field = Ktest*dense_density;
skeleton_field = Ktest(:,Jout)*(Scat_reference*Btest(Jin,:));
scattering_field = Ktest(:,Jout)*equivalent_density;
skeleton_field_error = relative_error(skeleton_field,exact_discrete_field);
scattering_field_error = relative_error(scattering_field,exact_discrete_field);
[exact_e,exact_h] = em3d.incoming_sources(zk,source,targets,'ehd');
analytic_field = [exact_e;exact_h];
finite_field_error = relative_error(reshape(scattering_field,6,[]),analytic_field);
dense_finite_error = relative_error(reshape(exact_discrete_field,6,[]),analytic_field);
fprintf('Held-out manufactured RHS error %.3e; discrete field compression %.3e\n', ...
    incoming_error,scattering_field_error)
assert(incoming_error<10*skeleton_tol)
assert(skeleton_field_error<200*skeleton_tol)
assert(scattering_field_error<100*factor_tol+200*skeleton_tol)

% Obtain amplitudes using the existing finite-target matrix only. Cancel the
% leading 1/R correction with two radii after removing exp(ikR)/R.
far_radii = [1e4,2e4];
far_scattering = complex(zeros(6*size(directions,2),2));
far_dense = far_scattering;
for j = 1:numel(far_radii)
    radius = far_radii(j);
    targets.r = radius*directions;
    Kfar = nrccie_matrix(S,targets,zk);
    scale = radius*exp(-1i*zk*radius);
    far_scattering(:,j) = scale*(Kfar(:,Jout)*equivalent_density);
    far_dense(:,j) = scale*(Kfar*dense_density);
end
farfield_scattering = reshape(2*far_scattering(:,2)-far_scattering(:,1),6,[]);
farfield_dense = reshape(2*far_dense(:,2)-far_dense(:,1),6,[]);

% Independent dipole asymptotics in the convention of incoming_sources:
% E_inf=-k^2 c [(I-dd^T)p+d x m], H_inf=k^2 c [(I-dd^T)m-d x p].
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
dense_farfield_e_error = relative_error(farfield_dense(1:3,:),electric_farfield);
dense_farfield_h_error = relative_error(farfield_dense(4:6,:),magnetic_farfield);
farfield_compression_error = relative_error(farfield_scattering,farfield_dense);
fprintf('Analytic farfield errors: E %.3e, H %.3e; dense E %.3e, H %.3e\n', ...
    farfield_e_error,farfield_h_error,dense_farfield_e_error,dense_farfield_h_error)
assert(farfield_compression_error<100*factor_tol+200*skeleton_tol)
% These looser bounds measure the coarse boundary discretization itself.
assert(max([finite_field_error,dense_finite_error,farfield_e_error, ...
    farfield_h_error,dense_farfield_e_error,dense_farfield_h_error])<3e-2)

% Confirm this test exercises low-rank elimination, not just root dense LU.
below_root = 0;
for j = 1:F.lvp(end-1)
    below_root = below_root+numel(F.factors(j).rd);
end
assert(below_root>0,'No low-rank elimination occurred before the root.')
factor_storage = whos('F');
fprintf('Surface: %d nodes, %d DOFs; FLAM levels %d, below-root eliminated %d\n', ...
    S.npts,ndof,F.nlvl,below_root)
fprintf('Factor: %.2f s, %.2f MiB (dense %.2f MiB); block error %.3e\n', ...
    factor_seconds,factor_storage.bytes/1024^2,16*ndof^2/1024^2,block_error)
fprintf('FMM/dense %.3e; rskelf_mv %.3e; rskelf_sv error %.3e, residual %.3e\n', ...
    fmm_dense_error,multiply_error,solve_error,solve_residual)
fprintf('Skeleton nodes: incoming %d, outgoing %d; Scat %d x %d, %.2f s\n', ...
    numel(Jin)/3,numel(Jout)/3,size(Scat,1),size(Scat,2),scattering_seconds)
fprintf('Scat/dense %.3e; held-out RHS %.3e, skeleton field %.3e, full field %.3e\n', ...
    scattering_error,incoming_error,skeleton_field_error,scattering_field_error)
fprintf('Analytic finite-distance field: scattering %.3e, dense %.3e\n', ...
    finite_field_error,dense_finite_error)
fprintf('Analytic farfield: E %.3e, H %.3e (dense E %.3e, H %.3e)\n', ...
    farfield_e_error,farfield_h_error,dense_farfield_e_error,dense_farfield_h_error)
fprintf('Farfield compression error %.3e; incoming training error %.3e\n', ...
    farfield_compression_error,incoming_training_error)
geometry = struct('semi_axes',[1.2,0.85,0.65],'patches_per_axis',[2,2,2], ...
    'center',[0;0;0],'order',4,'iptype',11);
settings = struct('zk',zk,'alpha',alpha,'jump',opts.jump, ...
    'eps_quad',eps_quad,'eps_fmm',fmm_opts.eps_fmm,'factor_tol',factor_tol, ...
    'skeleton_tol',skeleton_tol,'max_nodes',max_nodes,'max_patches',1, ...
    'batch_size',batch_size,'proxy_radius',4,'proxy_order',6, ...
    'training_positions',size(training_positions,2), ...
    'training_radius',0.25,'finite_target_radius',40,'far_radii',far_radii);
errors = struct('system_block',block_error,'fmm_dense',fmm_dense_error, ...
    'factor_multiply',multiply_error,'factor_solve',solve_error, ...
    'factor_residual',solve_residual,'scattering_dense',scattering_error, ...
    'incoming_heldout',incoming_error,'skeleton_field',skeleton_field_error, ...
    'scattering_field',scattering_field_error,'incoming_training',incoming_training_error, ...
    'finite_field',finite_field_error,'dense_finite_field',dense_finite_error, ...
    'farfield_e',farfield_e_error,'farfield_h',farfield_h_error, ...
    'dense_farfield_e',dense_farfield_e_error,'dense_farfield_h',dense_farfield_h_error, ...
    'farfield_compression',farfield_compression_error);
stats = struct('surface_nodes',S.npts,'factor_levels',F.nlvl, ...
    'eliminated_below_root',below_root,'factor_bytes',factor_storage.bytes, ...
    'dense_bytes',16*ndof^2,'factor_seconds',factor_seconds, ...
    'scattering_seconds',scattering_seconds,'incoming_nodes',numel(Jin)/3, ...
    'outgoing_nodes',numel(Jout)/3);
output_dir = fullfile(root,'diagnostics','nrccie_scattering');
if ~isfolder(output_dir), mkdir(output_dir); end
save(fullfile(output_dir,'scattering_results.mat'), ...
    'Scat','Jin','Jout','geometry','settings','errors','stats', ...
    'blocks_in','Pmerge_in','blocks_out','Pmerge_out', ...
    'source','training_positions','directions', ...
    'farfield_analytic','farfield_dense','farfield_scattering')
disp('NRCCIE_SCATTERING_MATRIX_PASS')


function rhs = manufactured_rhs(S,E,H,alpha)
normal_e = sum(S.n.*E,1);
tangent = -cross(S.n,H,1)+alpha*(S.n.*normal_e-E);
rhs = [sum(S.dru.*tangent,1);sum(S.drv.*tangent,1);-normal_e];
end


function err = relative_error(value,reference)
err = norm(value-reference,'fro')/max(norm(reference,'fro'),realmin);
end
