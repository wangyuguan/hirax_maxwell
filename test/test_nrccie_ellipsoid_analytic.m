% Server benchmark: 7,986-node ellipsoid using NRCCIE_APPLY.
% An interior dipole gives an analytic radiating exterior solution, while
% an exterior dipole gives an analytic regular interior solution. With the
% outward normal and this test's manufactured RHS and field representation,
% the interior combined trace uses jump=+0.5 and the exterior uses jump=-0.5.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir

%% Ellipsoid and matrix-free NRCCIE operator
semi_axes = [1.2,0.85,0.65];
patches_per_axis = [3,3,4];
surface_order = 10;
evaluation_order = 16;
zk = 0.9;
alpha = 1;
eps_quad = 1e-11;
eps_fmm = 1e-11;
eps_gmres = 1e-10;
maxit = 3;
restart = 100;

S = geometries.ellipsoid( ...
    semi_axes,patches_per_axis,[0;0;0],surface_order,11);
fprintf('Ellipsoid: %d patches, %d nodes, %d unknowns\n', ...
    S.npatches,S.npts,3*S.npts)
assert(abs(S.npts-8000)<=20, ...
    'This regression is intended to use approximately 8000 surface nodes.')
correction = nrccie_quad_corr_block(S,eps_quad,zk,[],256);
assert(sum(cellfun(@nnz,correction))>0)

% Oversampling is used only for independent smooth volume evaluation.
S_evaluation = oversample(S,evaluation_order);

cases = struct([]);
cases(1).name = 'exterior';
cases(1).source.r = [0.14;-0.09;0.07];
cases(1).source.edips = [1;0.25i;-0.3];
cases(1).source.hdips = [-0.2i;0.35;0.15i];
cases(1).probes = [ ...
     2.20,-2.00, 0.35, 1.55; ...
     0.20, 1.65,-1.80,-1.25; ...
     0.45,-0.75, 1.55, 1.30];

cases(2).name = 'interior';
cases(2).source.r = [2.6;-1.7;1.8];
cases(2).source.edips = [-0.4;0.7i;0.25];
cases(2).source.hdips = [0.3i;-0.2;0.5i];
cases(2).probes = [ ...
     0, 0.25,-0.35, 0.45; ...
     0, 0.15, 0.20,-0.25; ...
     0, 0.10,-0.15, 0.05];

results = struct([]);
for case_id = 1:numel(cases)
    test_case = cases(case_id);
    assert_probe_side(test_case.probes,semi_axes,test_case.name)

    [exact_boundary_e,exact_boundary_h] = em3d.incoming_sources( ...
        zk,test_case.source,S,'ehd');
    rhs_components = manufactured_rhs( ...
        S,exact_boundary_e,exact_boundary_h,alpha);
    rhs = rhs_components(:);

    options = struct('zk',zk,'alpha',alpha,'fmm',true, ...
        'eps_fmm',eps_fmm);
    if strcmp(test_case.name,'interior')
        options.jump = 0.5;
    else
        options.jump = -0.5;
    end
    matvec = @(density)nrccie_apply(S,density,correction,options);
    [~,explicit_matrix] = nrccie_apply( ...
        S,zeros(3*S.npts,1),correction,options);
    assert(isempty(explicit_matrix), ...
        'The iterative test unexpectedly formed a system matrix.')

    solve_timer = tic;
    [solution,flag,reported_residual,iterations,residual_history] = ...
        gmres(matvec,rhs,restart,eps_gmres,maxit);
    solve_time = toc(solve_timer);
    true_residual = norm(matvec(solution)-rhs)/norm(rhs);
    assert(flag==0 && true_residual<5*eps_gmres, ...
        '%s GMRES failed: flag %d, residual %.3e.', ...
        test_case.name,flag,true_residual)

    density_components = reshape(solution,3,S.npts);
    surface_current = S.dru.*density_components(1,:)+ ...
        S.drv.*density_components(2,:);
    cartesian_density = [surface_current;density_components(3,:)];
    evaluation_density = resample_values( ...
        S,cartesian_density,S_evaluation);
    [computed_e,computed_h] = direct_nrccie_fields( ...
        S_evaluation,evaluation_density,test_case.probes,zk);
    target.r = test_case.probes;
    [exact_e,exact_h] = em3d.incoming_sources( ...
        zk,test_case.source,target,'ehd');

    electric_error = relative_error(computed_e,exact_e);
    magnetic_error = relative_error(computed_h,exact_h);
    results(case_id).name = test_case.name;
    results(case_id).flag = flag;
    results(case_id).reported_residual = reported_residual;
    results(case_id).true_residual = true_residual;
    results(case_id).iterations = iterations;
    results(case_id).number_of_iterations = numel(residual_history)-1;
    results(case_id).solve_time = solve_time;
    results(case_id).electric_error = electric_error;
    results(case_id).magnetic_error = magnetic_error;

    fprintf(['%s: %d iterations, residual %.3e, ' ...
        'E error %.3e, H error %.3e, %.2f s\n'], ...
        test_case.name,results(case_id).number_of_iterations, ...
        true_residual,electric_error,magnetic_error,solve_time)
end

assert(max([results.electric_error])<5e-3, ...
    'Manufactured electric-field error is too large.')
assert(max([results.magnetic_error])<5e-3, ...
    'Manufactured magnetic-field error is too large.')
disp('NRCCIE_ELLIPSOID_ANALYTIC_PASS')


function rhs = manufactured_rhs(S,electric_field,magnetic_field,alpha)
normal_electric = sum(S.n.*electric_field,1);
nxnx_electric = S.n.*normal_electric-electric_field;
boundary_vector = -cross(S.n,magnetic_field,1)+alpha*nxnx_electric;
rhs = [sum(S.dru.*boundary_vector,1); ...
    sum(S.drv.*boundary_vector,1);-normal_electric];
end


function [electric_field,magnetic_field] = ...
    direct_nrccie_fields(S,density,targets,zk)
electric_field = complex(zeros(3,size(targets,2)));
magnetic_field = complex(zeros(3,size(targets,2)));
weights = S.wts(:).';
for target_id = 1:size(targets,2)
    displacement = targets(:,target_id)-S.r;
    distance = vecnorm(displacement,2,1);
    assert(all(distance>0),'A volume probe lies on a surface node.')
    green = exp(1i*zk*distance)./(4*pi*distance);
    radial = green.*(1i*zk*distance-1)./distance.^2;
    gradient = displacement.*radial;
    weighted_green = green.*weights;
    weighted_gradient = gradient.*weights;

    slp_current = sum(density(1:3,:).*weighted_green,2);
    gradient_charge = sum(weighted_gradient.*density(4,:),2);
    electric_field(:,target_id) = 1i*zk*slp_current-gradient_charge;
    magnetic_field(:,target_id) = [ ...
        sum(weighted_gradient(2,:).*density(3,:)- ...
            weighted_gradient(3,:).*density(2,:)); ...
        sum(weighted_gradient(3,:).*density(1,:)- ...
            weighted_gradient(1,:).*density(3,:)); ...
        sum(weighted_gradient(1,:).*density(2,:)- ...
            weighted_gradient(2,:).*density(1,:))];
end
end


function out = resample_values(S,values,T)
coefficients = vals2coefs(S,values);
out = complex(zeros(size(values,1),T.npts));
[kinds,~,groups] = unique([S.norders(:),S.iptype(:)],'rows');
for group = 1:size(kinds,1)
    patches = find(groups==group);
    norder = kinds(group,1);
    uv = T.uvs_targ(:,T.ixyzs(patches(1)):T.ixyzs(patches(1)+1)-1);
    if kinds(group,2)==1
        basis = koorn.pols(norder,uv);
    elseif kinds(group,2)==11
        basis = polytens.lege.pols(norder,uv);
    else
        basis = polytens.cheb.pols(norder,uv);
    end
    for patch = patches(:).'
        old = S.ixyzs(patch):S.ixyzs(patch+1)-1;
        new = T.ixyzs(patch):T.ixyzs(patch+1)-1;
        out(:,new) = coefficients(:,old)*basis;
    end
end
end


function assert_probe_side(points,semi_axes,side)
ellipsoid_level = sum((points./semi_axes(:)).^2,1);
if strcmp(side,'interior')
    assert(all(ellipsoid_level<0.8), ...
        'Interior probes are too close to or outside the ellipsoid.')
else
    assert(all(ellipsoid_level>1.4), ...
        'Exterior probes are too close to or inside the ellipsoid.')
end
end


function error = relative_error(computed,exact)
error = norm(computed(:)-exact(:))/norm(exact(:));
end
