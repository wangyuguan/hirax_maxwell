clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath('../FMM3D/matlab')

%% Configuration

sphere_radius = 1;
number_of_cube_subdivisions = 1;
surface_order = 8;
patch_type = 1;

zk = 1.1;
alpha = 1;
quadrature_tolerance = 1e-12;
fmm_tolerance = 1e-12;
gmres_tolerance = 1e-11;
maximum_iterations = 100;

required_residual = 5e-11;
required_exterior_field_error = 1e-8;

%% Unit-sphere surface and local tangent basis

S = geometries.sphere(sphere_radius,number_of_cube_subdivisions, ...
    [0;0;0],surface_order,patch_type);

surface_normal = S.n;
surface_ru = S.du./vecnorm(S.du,2,1);
surface_rv = cross(surface_normal,surface_ru,1);

fprintf('Independent Maxwell PEC NRCCIE sphere test\n')
fprintf('  patches: %d\n',S.npatches)
fprintf('  points: %d\n',S.npts)
fprintf('  surface order: %d\n',surface_order)

%% Add-subtract quadrature correction matrices

quadrature_timer = tic;
correction = struct();
[correction.slp,quadrature.slp] = em3d.slp.get_quad_corr_mat( ...
    S,quadrature_tolerance,zk);
[correction.dx,correction.dy,correction.dz,quadrature.sgrad] = ...
    em3d.sgrad.get_quad_corr_mat(S,quadrature_tolerance,zk);

quadrature.slp.kernel_order = 0;
oversampling_orders = get_oversampling_parameters( ...
    S,quadrature.slp,quadrature_tolerance);
[S_over,interpolation_to_over] = oversample(S,oversampling_orders);

correction = make_oversampled_corrections( ...
    S,S_over,interpolation_to_over,quadrature,zk);
quadrature_time = toc(quadrature_timer);

number_of_correction_entries = nnz(correction.slp)+ ...
    nnz(correction.dx)+nnz(correction.dy)+nnz(correction.dz);

fprintf('  correction entries: %d\n',number_of_correction_entries)
fprintf('  oversampled points: %d\n',S_over.npts)
fprintf('  correction construction: %.2f s\n',quadrature_time)

%% Analytic PEC multipole data and NRCCIE right-hand side

[electric_incident,magnetic_incident] = ...
    maxwell_multipole(S.r,zk,'regular');
[electric_outgoing_on_surface,magnetic_outgoing_on_surface] = ...
    maxwell_multipole(S.r,zk,'outgoing');

boundary_argument = zk*sphere_radius;
spherical_j1 = sin(boundary_argument)/boundary_argument^2- ...
    cos(boundary_argument)/boundary_argument;
spherical_y1 = -cos(boundary_argument)/boundary_argument^2- ...
    sin(boundary_argument)/boundary_argument;
scattering_coefficient = -spherical_j1/(spherical_j1+1i*spherical_y1);

electric_total_exact_on_surface = electric_incident+ ...
    scattering_coefficient*electric_outgoing_on_surface;
magnetic_total_exact_on_surface = magnetic_incident+ ...
    scattering_coefficient*magnetic_outgoing_on_surface;
surface_current_exact = cross(surface_normal, ...
    magnetic_total_exact_on_surface,1);
surface_charge_exact = sum(surface_normal.* ...
    electric_total_exact_on_surface,1);

normal_electric_incident = sum( ...
    surface_normal.*electric_incident,1);
normal_cross_magnetic_incident = cross( ...
    surface_normal,magnetic_incident,1);
normal_cross_normal_cross_electric_incident = ...
    surface_normal.*normal_electric_incident-electric_incident;

tangential_right_hand_side = normal_cross_magnetic_incident- ...
    alpha*normal_cross_normal_cross_electric_incident;

right_hand_side = complex(zeros(3,S.npts));
right_hand_side(1,:) = sum( ...
    surface_ru.*tangential_right_hand_side,1);
right_hand_side(2,:) = sum( ...
    surface_rv.*tangential_right_hand_side,1);
right_hand_side(3,:) = normal_electric_incident;
right_hand_side = right_hand_side(:);

%% Matrix-free NRCCIE solve using MATLAB GMRES and FMM3D

operator_data = struct();
operator_data.number_of_points = S.npts;
operator_data.target_points = S.r;
operator_data.source_points_over = S_over.r;
operator_data.surface_weights_over = S_over.wts(:).';
operator_data.interpolation_to_over = interpolation_to_over;
operator_data.surface_normal = surface_normal;
operator_data.surface_ru = surface_ru;
operator_data.surface_rv = surface_rv;
operator_data.zk = zk;
operator_data.alpha = alpha;
operator_data.fmm_tolerance = fmm_tolerance;
operator_data.correction = correction;

nrccie_operator = @(density) apply_nrccie(density,operator_data);

solve_timer = tic;
[solution,gmres_flag,~,~,residual_history] = gmres( ...
    nrccie_operator,right_hand_side,[],gmres_tolerance, ...
    maximum_iterations);
solve_time = toc(solve_timer);

number_of_iterations = length(residual_history)-1;
relative_residual = norm( ...
    nrccie_operator(solution)-right_hand_side)/norm(right_hand_side);

local_densities = reshape(solution,3,S.npts);
surface_current = surface_ru.*local_densities(1,:)+ ...
    surface_rv.*local_densities(2,:);
surface_charge = local_densities(3,:);
normal_current = max(abs(sum(surface_normal.*surface_current,1)));

local_densities_exact = [ ...
    sum(surface_ru.*surface_current_exact,1); ...
    sum(surface_rv.*surface_current_exact,1); ...
    surface_charge_exact];
exact_solution = local_densities_exact(:);
exact_solution_residual = norm( ...
    nrccie_operator(exact_solution)-right_hand_side)/ ...
    norm(right_hand_side);
density_difference = local_densities-local_densities_exact;
density_error = sqrt(sum(S.wts(:).'.*sum(abs(density_difference).^2,1)))/ ...
    sqrt(sum(S.wts(:).'.*sum(abs(local_densities_exact).^2,1)));

fprintf('  GMRES iterations: %d\n',number_of_iterations)
fprintf('  GMRES flag: %d\n',gmres_flag)
fprintf('  relative residual: %.3e\n',relative_residual)
fprintf('  maximum normal current: %.3e\n',normal_current)
fprintf('  surface-density error: %.3e\n',density_error)
fprintf('  analytic-density equation residual: %.3e\n', ...
    exact_solution_residual)
fprintf('  solve time: %.2f s\n',solve_time)

%% Analytic exterior scattered-field test

exterior_targets = [ ...
     1.45, -1.63,  0.71, -0.92,  1.31, -1.24; ...
     0.38,  0.62, -1.84,  1.57, -1.29, -0.88; ...
    -0.51,  0.79,  0.55, -0.74,  1.18, -1.46];

[electric_scattered,magnetic_scattered] = evaluate_nrccie_fields( ...
    surface_current,surface_charge,exterior_targets,operator_data);

target_info = struct();
target_info.r = exterior_targets;
[electric_outgoing,magnetic_outgoing] = ...
    maxwell_multipole(target_info.r,zk,'outgoing');

electric_exact = scattering_coefficient*electric_outgoing;
magnetic_exact = scattering_coefficient*magnetic_outgoing;

[electric_from_exact_density,magnetic_from_exact_density] = ...
    evaluate_nrccie_fields(surface_current_exact,surface_charge_exact, ...
    exterior_targets,operator_data);
exact_density_evaluation_error = norm( ...
    [electric_from_exact_density-electric_exact; ...
    magnetic_from_exact_density-magnetic_exact],'fro')/ ...
    norm([electric_exact;magnetic_exact],'fro');

exterior_field_error = norm( ...
    [electric_scattered-electric_exact; ...
    magnetic_scattered-magnetic_exact],'fro')/ ...
    norm([electric_exact;magnetic_exact],'fro');
exterior_digits = -log10(exterior_field_error);

fprintf('  exterior scattered-field error: %.3e\n', ...
    exterior_field_error)
fprintf('  exterior scattered-field digits: %.2f\n',exterior_digits)
fprintf('  exact-density exterior evaluation error: %.3e\n', ...
    exact_density_evaluation_error)

if gmres_flag == 0 && relative_residual < required_residual && ...
        exterior_field_error < required_exterior_field_error
    fprintf('  PASS\n')
else
    fprintf('  FAIL\n')
end


function correction = make_oversampled_corrections( ...
    S,S_over,interpolation_to_over,quadrature,zk)

number_of_points = S.npts;
[~,~,~,ixyzs] = extract_arrays(S);
[~,~,~,ixyzs_over] = extract_arrays(S_over);

accurate_slp = conv_rsc_to_spmat(S,quadrature.slp.row_ptr, ...
    quadrature.slp.col_ind,quadrature.slp.wnear);
accurate_dx = conv_rsc_to_spmat(S,quadrature.sgrad.row_ptr, ...
    quadrature.sgrad.col_ind,quadrature.sgrad.wnear(1,:).');
accurate_dy = conv_rsc_to_spmat(S,quadrature.sgrad.row_ptr, ...
    quadrature.sgrad.col_ind,quadrature.sgrad.wnear(2,:).');
accurate_dz = conv_rsc_to_spmat(S,quadrature.sgrad.row_ptr, ...
    quadrature.sgrad.col_ind,quadrature.sgrad.wnear(3,:).');

number_of_near_entries = nnz(accurate_slp);
row_indices = zeros(number_of_near_entries,1);
column_indices = zeros(number_of_near_entries,1);
smooth_slp_values = complex(zeros(number_of_near_entries,1));
smooth_dx_values = complex(zeros(number_of_near_entries,1));
smooth_dy_values = complex(zeros(number_of_near_entries,1));
smooth_dz_values = complex(zeros(number_of_near_entries,1));
entry_index = 0;

for target_index = 1:number_of_points
    interaction_indices = quadrature.slp.row_ptr(target_index): ...
        quadrature.slp.row_ptr(target_index+1)-1;
    for interaction_index = interaction_indices
        patch_index = quadrature.slp.col_ind(interaction_index);
        source_indices = ixyzs(patch_index):ixyzs(patch_index+1)-1;
        source_indices_over = ixyzs_over(patch_index): ...
            ixyzs_over(patch_index+1)-1;

        displacement = S.r(:,target_index)- ...
            S_over.r(:,source_indices_over);
        distance = vecnorm(displacement,2,1);
        use = distance > 1e-14;
        green_function = complex(zeros(size(distance)));
        radial_derivative = complex(zeros(size(distance)));
        green_function(use) = exp(1i*zk*distance(use))./ ...
            (4*pi*distance(use));
        radial_derivative(use) = (1i*zk*distance(use)-1).* ...
            exp(1i*zk*distance(use))./(4*pi*distance(use).^3);

        weights_over = reshape( ...
            S_over.wts(source_indices_over),[],1);
        weighted_interpolation = weights_over.* ...
            interpolation_to_over(source_indices_over,source_indices);

        new_entries = entry_index+(1:length(source_indices));
        row_indices(new_entries) = target_index;
        column_indices(new_entries) = source_indices;
        smooth_slp_values(new_entries) = ...
            green_function*weighted_interpolation;
        smooth_dx_values(new_entries) = ...
            (radial_derivative.*displacement(1,:))* ...
            weighted_interpolation;
        smooth_dy_values(new_entries) = ...
            (radial_derivative.*displacement(2,:))* ...
            weighted_interpolation;
        smooth_dz_values(new_entries) = ...
            (radial_derivative.*displacement(3,:))* ...
            weighted_interpolation;
        entry_index = entry_index+length(source_indices);
    end
end

row_indices = row_indices(1:entry_index);
column_indices = column_indices(1:entry_index);
smooth_slp = sparse(row_indices,column_indices, ...
    smooth_slp_values(1:entry_index),number_of_points,number_of_points);
smooth_dx = sparse(row_indices,column_indices, ...
    smooth_dx_values(1:entry_index),number_of_points,number_of_points);
smooth_dy = sparse(row_indices,column_indices, ...
    smooth_dy_values(1:entry_index),number_of_points,number_of_points);
smooth_dz = sparse(row_indices,column_indices, ...
    smooth_dz_values(1:entry_index),number_of_points,number_of_points);

correction = struct();
correction.slp = accurate_slp-smooth_slp;
correction.dx = accurate_dx-smooth_dx;
correction.dy = accurate_dy-smooth_dy;
correction.dz = accurate_dz-smooth_dz;
end


function [electric,magnetic] = maxwell_multipole(points,zk,kind)

radius = vecnorm(points,2,1);
argument = zk*radius;

spherical_j0 = sin(argument)./argument;
spherical_j1 = sin(argument)./argument.^2- ...
    cos(argument)./argument;

if strcmp(kind,'regular')
    radial_function = spherical_j1;
    radial_derivative = spherical_j0-2*spherical_j1./argument;
else
    spherical_y0 = -cos(argument)./argument;
    spherical_y1 = -cos(argument)./argument.^2- ...
        sin(argument)./argument;
    radial_function = spherical_j1+1i*spherical_y1;
    radial_derivative = spherical_j0+1i*spherical_y0- ...
        2*radial_function./argument;
end

x = points(1,:);
y = points(2,:);
z = points(3,:);
inverse_radius = 1./radius;
inverse_radius_squared = inverse_radius.^2;
cos_theta = z.*inverse_radius;
sin_theta_squared = (x.^2+y.^2).*inverse_radius_squared;

electric = radial_function.*[ ...
    -y.*inverse_radius; ...
     x.*inverse_radius; ...
     zeros(1,size(points,2))];

radial_vector = points.*inverse_radius;
sin_theta_theta_vector = [ ...
    z.*x.*inverse_radius_squared; ...
    z.*y.*inverse_radius_squared; ...
    -sin_theta_squared];

magnetic_n = 2*radial_function.*cos_theta./argument.* ...
    radial_vector-(radial_function+argument.*radial_derivative)./ ...
    argument.*sin_theta_theta_vector;
magnetic = -1i*magnetic_n;
end


function result = apply_nrccie(density,operator_data)

number_of_points = operator_data.number_of_points;
local_densities = reshape(density,3,number_of_points);

surface_current = operator_data.surface_ru.*local_densities(1,:)+ ...
    operator_data.surface_rv.*local_densities(2,:);
surface_charge = local_densities(3,:);
cartesian_densities = [surface_current;surface_charge];
cartesian_densities_over = (operator_data.interpolation_to_over* ...
    cartesian_densities.').';

fmm_sources = struct();
fmm_sources.sources = operator_data.source_points_over;
fmm_sources.nd = 4;
fmm_sources.charges = cartesian_densities_over.* ...
    operator_data.surface_weights_over;

fmm_output = hfmm3d(operator_data.fmm_tolerance, ...
    operator_data.zk,fmm_sources,0,operator_data.target_points,2);

potential = reshape(fmm_output.pottarg,4,number_of_points);
gradient = reshape(fmm_output.gradtarg,4,3,number_of_points);

potential = potential+(operator_data.correction.slp* ...
    cartesian_densities.').';
gradient_x = (operator_data.correction.dx*cartesian_densities.').';
gradient_y = (operator_data.correction.dy*cartesian_densities.').';
gradient_z = (operator_data.correction.dz*cartesian_densities.').';
gradient(:,1,:) = gradient(:,1,:)+reshape( ...
    gradient_x,4,1,number_of_points);
gradient(:,2,:) = gradient(:,2,:)+reshape( ...
    gradient_y,4,1,number_of_points);
gradient(:,3,:) = gradient(:,3,:)+reshape( ...
    gradient_z,4,1,number_of_points);

single_layer_current = potential(1:3,:);
single_layer_charge = potential(4,:);
gradient_charge = reshape(gradient(4,:,:),3,number_of_points);

curl_single_layer_current = complex(zeros(3,number_of_points));
curl_single_layer_current(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,number_of_points);
curl_single_layer_current(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,number_of_points);
curl_single_layer_current(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,number_of_points);

divergence_single_layer_current = reshape( ...
    gradient(1,1,:)+gradient(2,2,:)+gradient(3,3,:), ...
    1,number_of_points);

electric_scattered = 1i*operator_data.zk*single_layer_current- ...
    gradient_charge;
normal_cross_magnetic_scattered = cross( ...
    operator_data.surface_normal,curl_single_layer_current,1);
normal_electric_scattered = sum( ...
    operator_data.surface_normal.*electric_scattered,1);
normal_cross_normal_cross_electric_scattered = ...
    operator_data.surface_normal.*normal_electric_scattered- ...
    electric_scattered;

principal_value = complex(zeros(3,number_of_points));
first_equation = -normal_cross_magnetic_scattered+ ...
    operator_data.alpha* ...
    normal_cross_normal_cross_electric_scattered;
principal_value(1,:) = sum( ...
    operator_data.surface_ru.*first_equation,1);
principal_value(2,:) = sum( ...
    operator_data.surface_rv.*first_equation,1);
principal_value(3,:) = -normal_electric_scattered+ ...
    operator_data.alpha*(divergence_single_layer_current- ...
    1i*operator_data.zk*single_layer_charge);

result = (0.5*local_densities+principal_value);
result = result(:);
end


function [electric,magnetic] = evaluate_nrccie_fields( ...
    surface_current,surface_charge,target_points,operator_data)

number_of_targets = size(target_points,2);
cartesian_densities = [surface_current;surface_charge];
cartesian_densities_over = (operator_data.interpolation_to_over* ...
    cartesian_densities.').';

fmm_sources = struct();
fmm_sources.sources = operator_data.source_points_over;
fmm_sources.nd = 4;
fmm_sources.charges = cartesian_densities_over.* ...
    operator_data.surface_weights_over;

fmm_output = hfmm3d(operator_data.fmm_tolerance, ...
    operator_data.zk,fmm_sources,0,target_points,2);

potential = reshape(fmm_output.pottarg,4,number_of_targets);
gradient = reshape(fmm_output.gradtarg,4,3,number_of_targets);

electric = 1i*operator_data.zk*potential(1:3,:)- ...
    reshape(gradient(4,:,:),3,number_of_targets);

magnetic = complex(zeros(3,number_of_targets));
magnetic(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,number_of_targets);
magnetic(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,number_of_targets);
magnetic(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,number_of_targets);
end
