% Test complete NRCCIE surface-to-surface operator blocks analytically.

clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath('../FMM3D/matlab')

%% Sphere geometry and NRCCIE parameters

radius = 1;
surface_gap = radius/30;
center_offset = radius+surface_gap/2;
source_center = [-center_offset; center_offset;0];
target_center = [ center_offset; center_offset;0];
lower_left_center = [-center_offset;-center_offset;0];
lower_right_center = [ center_offset;-center_offset;0];

norder = 10;
sphere_subdivisions = 2;
sphere_iptype = 1;

wavelength = 10*radius;
zk = 2*pi/wavelength;
alpha = 1;
eps_quad = 1e-10;
eps_fmm = 1e-12;
maximum_spherical_harmonic_degree = 2;

source_surface = geometries.sphere( ...
    radius,sphere_subdivisions,source_center,norder,sphere_iptype);
close_target_surface = geometries.sphere( ...
    radius,sphere_subdivisions,target_center,norder,sphere_iptype);
lower_left_surface = geometries.sphere( ...
    radius,sphere_subdivisions,lower_left_center,norder,sphere_iptype);
lower_right_surface = geometries.sphere( ...
    radius,sphere_subdivisions,lower_right_center,norder,sphere_iptype);

target_surfaces = {source_surface,close_target_surface};
mapping_names = {'self 1 <- 1','close 2 <- 1'};
system_surfaces = {source_surface,close_target_surface, ...
    lower_left_surface,lower_right_surface};
system_centers = {source_center,target_center, ...
    lower_left_center,lower_right_center};
system_surface = merge([source_surface,close_target_surface, ...
    lower_left_surface,lower_right_surface]);

source_normal = source_surface.n;
source_ru = source_surface.du./vecnorm(source_surface.du,2,1);
source_rv = cross(source_normal,source_ru,1);

%% Charge, poloidal-current, and toroidal-current test modes

number_of_tests = 3*(maximum_spherical_harmonic_degree+1)^2-2;
test_type = cell(1,number_of_tests);
test_degree = zeros(1,number_of_tests);
test_order = zeros(1,number_of_tests);
test_count = 0;

for ell = 0:maximum_spherical_harmonic_degree
    for emm = -ell:ell
        test_count = test_count+1;
        test_type{test_count} = 'charge';
        test_degree(test_count) = ell;
        test_order(test_count) = emm;

        if ell>0
            test_count = test_count+1;
            test_type{test_count} = 'poloidal';
            test_degree(test_count) = ell;
            test_order(test_count) = emm;

            test_count = test_count+1;
            test_type{test_count} = 'toroidal';
            test_degree(test_count) = ell;
            test_order(test_count) = emm;
        end
    end
end

source_components = complex(zeros(3,source_surface.npts,number_of_tests));
source_cartesian_densities = ...
    complex(zeros(4,source_surface.npts,number_of_tests));
source_relative_position = source_surface.r-source_center;

for test = 1:number_of_tests
    [Y,surface_gradient_Y] = spherical_harmonic( ...
        test_degree(test),test_order(test),source_relative_position);
    surface_current = complex(zeros(3,source_surface.npts));
    surface_charge = complex(zeros(1,source_surface.npts));

    if strcmp(test_type{test},'charge')
        surface_charge = Y;
    elseif strcmp(test_type{test},'poloidal')
        surface_current = surface_gradient_Y/radius;
    else
        surface_current = cross( ...
            source_normal,surface_gradient_Y/radius,1);
    end

    current_u = sum(source_ru.*surface_current,1);
    current_v = sum(source_rv.*surface_current,1);
    surface_current = source_ru.*current_u+source_rv.*current_v;

    source_components(:,:,test) = [current_u;current_v;surface_charge];
    source_cartesian_densities(:,:,test) = ...
        [surface_current;surface_charge];
end

%% Independent spherical-harmonic expansion of the vector densities

maximum_expansion_degree = maximum_spherical_harmonic_degree+1;
number_of_expansion_modes = (maximum_expansion_degree+1)^2;
expansion_degree = zeros(1,number_of_expansion_modes);
expansion_order = zeros(1,number_of_expansion_modes);
expansion_count = 0;

for ell = 0:maximum_expansion_degree
    for emm = -ell:ell
        expansion_count = expansion_count+1;
        expansion_degree(expansion_count) = ell;
        expansion_order(expansion_count) = emm;
    end
end

test_expansion_mode_index = zeros(1,number_of_tests);
for test = 1:number_of_tests
    test_expansion_mode_index(test) = find( ...
        expansion_degree==test_degree(test) & ...
        expansion_order==test_order(test),1);
end

number_of_transform_theta_nodes = 24;
number_of_transform_phi_nodes = 48;
[cos_theta_nodes,cos_theta_weights] = ...
    gauss_legendre(number_of_transform_theta_nodes);
phi_nodes = 2*pi*(0:number_of_transform_phi_nodes-1)/ ...
    number_of_transform_phi_nodes;
[cos_theta_grid,phi_grid] = ndgrid(cos_theta_nodes,phi_nodes);
sin_theta_grid = sqrt(1-cos_theta_grid.^2);
transform_directions = [ ...
    reshape(sin_theta_grid.*cos(phi_grid),1,[]); ...
    reshape(sin_theta_grid.*sin(phi_grid),1,[]); ...
    reshape(cos_theta_grid,1,[])];
transform_weights = reshape( ...
    repmat(cos_theta_weights(:),1,number_of_transform_phi_nodes)* ...
    (2*pi/number_of_transform_phi_nodes),1,[]);

transform_harmonics = complex(zeros( ...
    number_of_expansion_modes,size(transform_directions,2)));
source_harmonics = complex(zeros( ...
    number_of_expansion_modes,source_surface.npts));

for mode = 1:number_of_expansion_modes
    transform_harmonics(mode,:) = spherical_harmonic( ...
        expansion_degree(mode),expansion_order(mode),transform_directions);
    source_harmonics(mode,:) = spherical_harmonic( ...
        expansion_degree(mode),expansion_order(mode), ...
        source_relative_position);
end

spectral_coefficients = complex(zeros( ...
    4,number_of_expansion_modes,number_of_tests));
spectral_reconstruction_error = zeros(1,number_of_tests);

for test = 1:number_of_tests
    [Y,surface_gradient_Y] = spherical_harmonic( ...
        test_degree(test),test_order(test),transform_directions);
    current_on_transform_grid = complex(zeros( ...
        3,size(transform_directions,2)));
    charge_on_transform_grid = complex(zeros( ...
        1,size(transform_directions,2)));

    if strcmp(test_type{test},'charge')
        charge_on_transform_grid = Y;
    elseif strcmp(test_type{test},'poloidal')
        current_on_transform_grid = surface_gradient_Y/radius;
    else
        current_on_transform_grid = cross( ...
            transform_directions,surface_gradient_Y/radius,1);
    end

    density_on_transform_grid = ...
        [current_on_transform_grid;charge_on_transform_grid];
    coefficients = (density_on_transform_grid.*transform_weights)* ...
        conj(transform_harmonics).';
    spectral_coefficients(:,:,test) = coefficients;

    reconstructed_density = coefficients*source_harmonics;
    original_density = source_cartesian_densities(:,:,test);
    reconstruction_difference = reconstructed_density-original_density;
    source_weights = source_surface.wts(:).';
    spectral_reconstruction_error(test) = sqrt(sum(source_weights.*sum( ...
        abs(reconstruction_difference).^2,1)))/sqrt(sum( ...
        source_weights.*sum(abs(original_density).^2,1)));
end

%% Complete NRCCIE self and close cross-surface mappings

smooth_operator_error = zeros(2,number_of_tests);
corrected_operator_error = zeros(2,number_of_tests);
corrected_operator_absolute_error = zeros(2,number_of_tests);
operator_reference_norm = zeros(2,number_of_tests);
correction_nonzeros = zeros(2,4);

for mapping = 1:2
    target_surface = target_surfaces{mapping};
    self_interaction = mapping==1;
    number_of_targets = target_surface.npts;
    target_normal = target_surface.n;
    target_ru = target_surface.du./vecnorm(target_surface.du,2,1);
    target_rv = cross(target_normal,target_ru,1);
    target_weights = target_surface.wts(:).';

    if self_interaction
        Cslp = em3d.slp.get_quad_corr_mat( ...
            source_surface,eps_quad,zk);
        [Cx,Cy,Cz] = em3d.sgrad.get_quad_corr_mat( ...
            source_surface,eps_quad,zk);
    else
        Cslp = em3d.slp.get_quad_corr_mat( ...
            source_surface,eps_quad,zk,target_surface);
        [Cx,Cy,Cz] = em3d.sgrad.get_quad_corr_mat( ...
            source_surface,eps_quad,zk,target_surface);
    end

    correction_nonzeros(mapping,:) = [nnz(Cslp),nnz(Cx),nnz(Cy),nnz(Cz)];

    analytic_scalar_potential = complex(zeros( ...
        number_of_expansion_modes,number_of_targets));
    analytic_scalar_gradient = complex(zeros( ...
        number_of_expansion_modes,3,number_of_targets));

    for mode = 1:number_of_expansion_modes
        [potential_reference,gradient_reference] = ...
            spherical_harmonic_slp_reference( ...
            expansion_degree(mode),expansion_order(mode), ...
            target_surface.r,source_center,radius,zk,self_interaction);
        analytic_scalar_potential(mode,:) = potential_reference;
        analytic_scalar_gradient(mode,:,:) = reshape( ...
            gradient_reference,1,3,number_of_targets);
    end

    for test = 1:number_of_tests
        cartesian_density = source_cartesian_densities(:,:,test);
        source = struct();
        source.sources = source_surface.r;
        source.nd = 4;
        source.charges = cartesian_density.*source_surface.wts(:).';

        if self_interaction
            fmm_output = hfmm3d(eps_fmm,zk,source,2);
            potential_smooth = reshape(fmm_output.pot,4,number_of_targets);
            gradient_smooth = reshape(fmm_output.grad,4,3,number_of_targets);
            identity_components = source_components(:,:,test);
        else
            fmm_output = hfmm3d( ...
                eps_fmm,zk,source,0,target_surface.r,2);
            potential_smooth = reshape( ...
                fmm_output.pottarg,4,number_of_targets);
            gradient_smooth = reshape( ...
                fmm_output.gradtarg,4,3,number_of_targets);
            identity_components = complex(zeros(3,number_of_targets));
        end

        potential_corrected = potential_smooth+ ...
            (Cslp*cartesian_density.').';
        gradient_corrected = gradient_smooth;
        gradient_corrected(:,1,:) = gradient_corrected(:,1,:)+reshape( ...
            (Cx*cartesian_density.').',4,1,number_of_targets);
        gradient_corrected(:,2,:) = gradient_corrected(:,2,:)+reshape( ...
            (Cy*cartesian_density.').',4,1,number_of_targets);
        gradient_corrected(:,3,:) = gradient_corrected(:,3,:)+reshape( ...
            (Cz*cartesian_density.').',4,1,number_of_targets);

        coefficients = spectral_coefficients(:,:,test);
        potential_reference = coefficients*analytic_scalar_potential;
        gradient_reference = complex(zeros(4,3,number_of_targets));
        for direction = 1:3
            scalar_gradient_direction = reshape( ...
                analytic_scalar_gradient(:,direction,:), ...
                number_of_expansion_modes,number_of_targets);
            gradient_reference(:,direction,:) = reshape( ...
                coefficients*scalar_gradient_direction, ...
                4,1,number_of_targets);
        end

        operator_smooth = assemble_nrccie_numerical( ...
            potential_smooth,gradient_smooth,target_normal,target_ru, ...
            target_rv,zk,alpha,identity_components,self_interaction);
        operator_corrected = assemble_nrccie_numerical( ...
            potential_corrected,gradient_corrected,target_normal,target_ru, ...
            target_rv,zk,alpha,identity_components,self_interaction);
        scalar_harmonic_potential = analytic_scalar_potential( ...
            test_expansion_mode_index(test),:);
        operator_reference = assemble_nrccie_family_reference( ...
            potential_reference,gradient_reference, ...
            scalar_harmonic_potential,target_normal,target_ru,target_rv, ...
            zk,alpha,identity_components,self_interaction, ...
            test_type{test},test_degree(test),radius);

        reference_norm = sqrt(sum(target_weights.*sum( ...
            abs(operator_reference).^2,1)));
        operator_reference_norm(mapping,test) = reference_norm;
        smooth_operator_error(mapping,test) = sqrt(sum(target_weights.*sum( ...
            abs(operator_smooth-operator_reference).^2,1)))/reference_norm;
        corrected_difference_norm = sqrt(sum( ...
            target_weights.*sum(abs( ...
            operator_corrected-operator_reference).^2,1)));
        corrected_operator_absolute_error(mapping,test) = ...
            corrected_difference_norm;
        corrected_operator_error(mapping,test) = ...
            corrected_difference_norm/reference_norm;
    end

    clear Cslp Cx Cy Cz
end

%% Complete four-surface NRCCIE system

number_of_surfaces = 4;
nodes_per_surface = source_surface.npts;
surface_amplitude_cases = eye(number_of_surfaces);
number_of_amplitude_cases = size(surface_amplitude_cases,2);

system_normal = system_surface.n;
system_ru = system_surface.du./vecnorm(system_surface.du,2,1);
system_rv = cross(system_normal,system_ru,1);
system_weights = system_surface.wts(:).';

Cslp_system = em3d.slp.get_quad_corr_mat( ...
    system_surface,eps_quad,zk);
[Cx_system,Cy_system,Cz_system] = em3d.sgrad.get_quad_corr_mat( ...
    system_surface,eps_quad,zk);

system_correction_blocks = zeros(number_of_surfaces);
for target_surface_number = 1:number_of_surfaces
    target_ids = (target_surface_number-1)*nodes_per_surface+ ...
        (1:nodes_per_surface);
    for source_surface_number = 1:number_of_surfaces
        source_ids = (source_surface_number-1)*nodes_per_surface+ ...
            (1:nodes_per_surface);
        system_correction_blocks(target_surface_number,source_surface_number) = ...
            nnz(Cslp_system(target_ids,source_ids));
    end
end

system_scalar_potential = cell(number_of_surfaces,number_of_surfaces);
system_scalar_gradient = cell(number_of_surfaces,number_of_surfaces);

for target_surface_number = 1:number_of_surfaces
    target_surface = system_surfaces{target_surface_number};
    for source_surface_number = 1:number_of_surfaces
        self_interaction = target_surface_number==source_surface_number;
        source_center_use = system_centers{source_surface_number};
        scalar_potential = complex(zeros( ...
            number_of_expansion_modes,nodes_per_surface));
        scalar_gradient = complex(zeros( ...
            number_of_expansion_modes,3,nodes_per_surface));

        for mode = 1:number_of_expansion_modes
            [potential_reference,gradient_reference] = ...
                spherical_harmonic_slp_reference( ...
                expansion_degree(mode),expansion_order(mode), ...
                target_surface.r,source_center_use,radius,zk, ...
                self_interaction);
            scalar_potential(mode,:) = potential_reference;
            scalar_gradient(mode,:,:) = reshape( ...
                gradient_reference,1,3,nodes_per_surface);
        end

        system_scalar_potential{target_surface_number,source_surface_number} = ...
            scalar_potential;
        system_scalar_gradient{target_surface_number,source_surface_number} = ...
            scalar_gradient;
    end
end

system_smooth_error = zeros(number_of_amplitude_cases,number_of_tests);
system_corrected_error = zeros(number_of_amplitude_cases,number_of_tests);
system_corrected_absolute_error = zeros( ...
    number_of_amplitude_cases,number_of_tests);
system_reference_norm = zeros(number_of_amplitude_cases,number_of_tests);
system_corrected_block_error = zeros( ...
    number_of_surfaces,number_of_surfaces,number_of_tests);

for amplitude_case = 1:number_of_amplitude_cases
  surface_amplitudes = surface_amplitude_cases(:,amplitude_case);
  for test = 1:number_of_tests
    system_density = complex(zeros(4,system_surface.npts));
    system_identity_components = complex(zeros(3,system_surface.npts));

    for source_surface_number = 1:number_of_surfaces
        source_ids = (source_surface_number-1)*nodes_per_surface+ ...
            (1:nodes_per_surface);
        system_density(:,source_ids) = surface_amplitudes( ...
            source_surface_number)*source_cartesian_densities(:,:,test);
        system_identity_components(:,source_ids) = surface_amplitudes( ...
            source_surface_number)*source_components(:,:,test);
    end

    source = struct();
    source.sources = system_surface.r;
    source.nd = 4;
    source.charges = system_density.*system_surface.wts(:).';
    fmm_output = hfmm3d(eps_fmm,zk,source,2);
    potential_smooth = reshape(fmm_output.pot,4,system_surface.npts);
    gradient_smooth = reshape(fmm_output.grad,4,3,system_surface.npts);

    potential_corrected = potential_smooth+ ...
        (Cslp_system*system_density.').';
    gradient_corrected = gradient_smooth;
    gradient_corrected(:,1,:) = gradient_corrected(:,1,:)+reshape( ...
        (Cx_system*system_density.').',4,1,system_surface.npts);
    gradient_corrected(:,2,:) = gradient_corrected(:,2,:)+reshape( ...
        (Cy_system*system_density.').',4,1,system_surface.npts);
    gradient_corrected(:,3,:) = gradient_corrected(:,3,:)+reshape( ...
        (Cz_system*system_density.').',4,1,system_surface.npts);

    potential_reference = complex(zeros(4,system_surface.npts));
    gradient_reference = complex(zeros(4,3,system_surface.npts));
    scalar_harmonic_potential_reference = complex( ...
        zeros(1,system_surface.npts));

    for target_surface_number = 1:number_of_surfaces
        target_ids = (target_surface_number-1)*nodes_per_surface+ ...
            (1:nodes_per_surface);
        potential_target = complex(zeros(4,nodes_per_surface));
        gradient_target = complex(zeros(4,3,nodes_per_surface));
        scalar_harmonic_potential_target = complex( ...
            zeros(1,nodes_per_surface));

        for source_surface_number = 1:number_of_surfaces
            coefficients = surface_amplitudes(source_surface_number)* ...
                spectral_coefficients(:,:,test);
            potential_target = potential_target+coefficients* ...
                system_scalar_potential{ ...
                target_surface_number,source_surface_number};
            scalar_harmonic_potential_target = ...
                scalar_harmonic_potential_target+ ...
                surface_amplitudes(source_surface_number)* ...
                system_scalar_potential{ ...
                target_surface_number,source_surface_number}( ...
                test_expansion_mode_index(test),:);

            scalar_gradient = system_scalar_gradient{ ...
                target_surface_number,source_surface_number};
            for direction = 1:3
                scalar_gradient_direction = reshape( ...
                    scalar_gradient(:,direction,:), ...
                    number_of_expansion_modes,nodes_per_surface);
                gradient_target(:,direction,:) = ...
                    gradient_target(:,direction,:)+reshape( ...
                    coefficients*scalar_gradient_direction, ...
                    4,1,nodes_per_surface);
            end
        end

        potential_reference(:,target_ids) = potential_target;
        gradient_reference(:,:,target_ids) = gradient_target;
        scalar_harmonic_potential_reference(:,target_ids) = ...
            scalar_harmonic_potential_target;
    end

    operator_smooth = assemble_nrccie_numerical( ...
        potential_smooth,gradient_smooth,system_normal,system_ru, ...
        system_rv,zk,alpha,system_identity_components,true);
    operator_corrected = assemble_nrccie_numerical( ...
        potential_corrected,gradient_corrected,system_normal,system_ru, ...
        system_rv,zk,alpha,system_identity_components,true);
    operator_reference = assemble_nrccie_family_reference( ...
        potential_reference,gradient_reference, ...
        scalar_harmonic_potential_reference,system_normal,system_ru, ...
        system_rv,zk,alpha,system_identity_components,true, ...
        test_type{test},test_degree(test),radius);

    for target_surface_number = 1:number_of_surfaces
        target_ids = (target_surface_number-1)*nodes_per_surface+ ...
            (1:nodes_per_surface);
        target_weights = system_weights(target_ids);
        block_reference = operator_reference(:,target_ids);
        block_difference = operator_corrected(:,target_ids)- ...
            block_reference;
        block_reference_norm = sqrt(sum(target_weights.*sum( ...
            abs(block_reference).^2,1)));
        system_corrected_block_error( ...
            target_surface_number,amplitude_case,test) = sqrt(sum( ...
            target_weights.*sum(abs(block_difference).^2,1)))/ ...
            block_reference_norm;
    end

    reference_norm = sqrt(sum(system_weights.*sum( ...
        abs(operator_reference).^2,1)));
    system_reference_norm(amplitude_case,test) = reference_norm;
    system_smooth_error(amplitude_case,test) = sqrt(sum( ...
        system_weights.*sum( ...
        abs(operator_smooth-operator_reference).^2,1)))/reference_norm;
    corrected_difference_norm = sqrt(sum( ...
        system_weights.*sum(abs( ...
        operator_corrected-operator_reference).^2,1)));
    system_corrected_absolute_error(amplitude_case,test) = ...
        corrected_difference_norm;
    system_corrected_error(amplitude_case,test) = ...
        corrected_difference_norm/reference_norm;
  end
end

%% Results

fprintf('Analytic NRCCIE surface-operator test\n')
fprintf('  radius: %.8g, close-surface gap: %.8g = R/30\n', ...
    radius,surface_gap)
fprintf('  surface order: %d, nodes per sphere: %d\n', ...
    norder,source_surface.npts)
fprintf('  input harmonics: 0 <= ell <= %d\n', ...
    maximum_spherical_harmonic_degree)
fprintf('  maximum independent spectral reconstruction error: %.3e\n\n', ...
    max(spectral_reconstruction_error))

fprintf('mapping          smooth NRCCIE    corrected NRCCIE    nnz S / dx / dy / dz\n')
for mapping = 1:2
    fprintf('%-13s    %11.3e       %11.3e       %d / %d / %d / %d\n', ...
        mapping_names{mapping},max(smooth_operator_error(mapping,:)), ...
        max(corrected_operator_error(mapping,:)), ...
        correction_nonzeros(mapping,1),correction_nonzeros(mapping,2), ...
        correction_nonzeros(mapping,3),correction_nonzeros(mapping,4))
end

fprintf('\nFour-surface global NRCCIE system\n')
fprintf(['  four independent source-surface inputs sample all surface ' ...
    'block columns\n'])
fprintf('  maximum smooth error: %.3e\n',max(system_smooth_error,[],'all'))
fprintf('  maximum corrected error: %.3e\n', ...
    max(system_corrected_error,[],'all'))
for amplitude_case = 1:number_of_amplitude_cases
    fprintf('  source sphere %d maximum corrected error: %.3e\n', ...
        amplitude_case,max(system_corrected_error(amplitude_case,:)))
end
fprintf('  S correction blocks; rows target, columns source:\n')
disp(system_correction_blocks)
fprintf('  maximum corrected relative error in each individual block:\n')
maximum_system_block_error = max(system_corrected_block_error,[],3);
for target_surface_number = 1:number_of_surfaces
    for source_surface_number = 1:number_of_surfaces
        fprintf(' %11.3e',maximum_system_block_error( ...
            target_surface_number,source_surface_number))
    end
    fprintf('\n')
end

family_names = {'charge','poloidal','toroidal'};
fprintf('Corrected L2 results by density family\n')
fprintf(['family          self relative     close relative    system relative  ' ...
    ' system absolute   minimum system norm\n'])
for family = 1:numel(family_names)
    family_ids = strcmp(test_type,family_names{family});
    fprintf('%-10s      %11.3e       %11.3e       %11.3e    %11.3e       %11.3e\n', ...
        family_names{family}, ...
        max(corrected_operator_error(1,family_ids)), ...
        max(corrected_operator_error(2,family_ids)), ...
        max(system_corrected_error(:,family_ids),[],'all'), ...
        max(system_corrected_absolute_error(:,family_ids),[],'all'), ...
        min(system_reference_norm(:,family_ids),[],'all'))
end

fprintf('\nCorrected full-NRCCIE relative L2 errors by input mode\n')
fprintf(['type          ell    m       self 1 <- 1      close 2 <- 1    ' ...
    'four-surface system\n'])
for test = 1:number_of_tests
    fprintf('%-10s    %2d   %3d       %11.3e       %11.3e       %11.3e\n', ...
        test_type{test},test_degree(test),test_order(test), ...
        corrected_operator_error(1,test),corrected_operator_error(2,test), ...
        max(system_corrected_error(:,test)))
end


function values = assemble_nrccie_numerical( ...
    potential,gradient,normal,ru,rv,zk,alpha,identity_components, ...
    self_interaction)

number_of_targets = size(potential,2);
slp_current = potential(1:3,:);
slp_charge = potential(4,:);
gradient_slp_charge = reshape(gradient(4,:,:),3,number_of_targets);

curl_slp_current = complex(zeros(3,number_of_targets));
curl_slp_current(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,number_of_targets);
curl_slp_current(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,number_of_targets);
curl_slp_current(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,number_of_targets);
divergence_slp_current = reshape( ...
    gradient(1,1,:)+gradient(2,2,:)+gradient(3,3,:),1,number_of_targets);

electric_field = 1i*zk*slp_current-gradient_slp_charge;
normal_electric_field = sum(normal.*electric_field,1);
nxnxe = normal.*normal_electric_field-electric_field;
tangent_equation = -cross(normal,curl_slp_current,1)+alpha*nxnxe;

values = complex(zeros(3,number_of_targets));
values(1,:) = sum(ru.*tangent_equation,1);
values(2,:) = sum(rv.*tangent_equation,1);
values(3,:) = -normal_electric_field+alpha*( ...
    divergence_slp_current-1i*zk*slp_charge);

if self_interaction
    values = values+0.5*identity_components;
end
end


function values = assemble_nrccie_family_reference( ...
    potential,gradient,scalar_harmonic_potential,normal,ru,rv,zk, ...
    alpha,identity_components,self_interaction,test_type,ell,radius)

number_of_targets = size(potential,2);

if strcmp(test_type,'charge')
    gradient_scalar_potential = reshape( ...
        gradient(4,:,:),3,number_of_targets);
    normal_gradient_scalar_potential = sum( ...
        normal.*gradient_scalar_potential,1);
    tangent_equation = alpha*(gradient_scalar_potential- ...
        normal.*normal_gradient_scalar_potential);
    scalar_equation = normal_gradient_scalar_potential- ...
        1i*alpha*zk*scalar_harmonic_potential;
else
    vector_potential = potential(1:3,:);
    normal_vector_potential = sum(normal.*vector_potential,1);
    nxnx_vector_potential = normal.*normal_vector_potential- ...
        vector_potential;
    curl_vector_potential = complex(zeros(3,number_of_targets));
    curl_vector_potential(1,:) = reshape( ...
        gradient(3,2,:)-gradient(2,3,:),1,number_of_targets);
    curl_vector_potential(2,:) = reshape( ...
        gradient(1,3,:)-gradient(3,1,:),1,number_of_targets);
    curl_vector_potential(3,:) = reshape( ...
        gradient(2,1,:)-gradient(1,2,:),1,number_of_targets);
    tangent_equation = -cross(normal,curl_vector_potential,1)+ ...
        1i*alpha*zk*nxnx_vector_potential;
    scalar_equation = -1i*zk*normal_vector_potential;

    if strcmp(test_type,'poloidal')
        scalar_equation = scalar_equation- ...
            alpha*ell*(ell+1)/radius^2* ...
            scalar_harmonic_potential;
    end
end

values = complex(zeros(3,number_of_targets));
values(1,:) = sum(ru.*tangent_equation,1);
values(2,:) = sum(rv.*tangent_equation,1);
values(3,:) = scalar_equation;

if self_interaction
    values = values+0.5*identity_components;
end
end


function [potential,gradient] = spherical_harmonic_slp_reference( ...
    ell,emm,target_points,source_center,radius,zk,self_interaction)

relative_position = target_points-source_center;
radial_distance = vecnorm(relative_position,2,1);
radial_direction = relative_position./radial_distance;
[Y,surface_gradient_Y] = spherical_harmonic( ...
    ell,emm,relative_position);

zradius = zk*radius;
jell = sqrt(pi/(2*zradius))*besselj(ell+0.5,zradius);
jell_next = sqrt(pi/(2*zradius))*besselj(ell+1.5,zradius);
hell = sqrt(pi/(2*zradius))*besselh(ell+0.5,1,zradius);
hell_next = sqrt(pi/(2*zradius))*besselh(ell+1.5,1,zradius);
jell_derivative = ell/zradius*jell-jell_next;
hell_derivative = ell/zradius*hell-hell_next;

if self_interaction
    radial_factor = 1i*zk*radius^2*jell*hell;
    radial_derivative = 0.5i*zk^2*radius^2*( ...
        jell*hell_derivative+jell_derivative*hell);
    radial_distance(:) = radius;
else
    ztarget = zk*radial_distance;
    htarget = sqrt(pi./(2*ztarget)).* ...
        besselh(ell+0.5,1,ztarget);
    htarget_next = sqrt(pi./(2*ztarget)).* ...
        besselh(ell+1.5,1,ztarget);
    htarget_derivative = ell./ztarget.*htarget-htarget_next;
    radial_factor = 1i*zk*radius^2*jell.*htarget;
    radial_derivative = 1i*zk^2*radius^2*jell.*htarget_derivative;
end

potential = radial_factor.*Y;
gradient = radial_direction.*(radial_derivative.*Y)+ ...
    surface_gradient_Y.*(radial_factor./radial_distance);
end


function [Y,surface_gradient_Y] = spherical_harmonic( ...
    ell,emm,relative_position)

radial_distance = vecnorm(relative_position,2,1);
direction = relative_position./radial_distance;
cos_theta = direction(3,:);
sin_theta = sqrt(direction(1,:).^2+direction(2,:).^2);
phi = atan2(direction(2,:),direction(1,:));
absolute_order = abs(emm);

associated_values = legendre(ell,cos_theta,'unnorm');
associated_value = reshape( ...
    associated_values(absolute_order+1,:),1,[]);

if ell==0
    associated_theta_derivative = zeros(size(cos_theta));
else
    if absolute_order<=ell-1
        previous_values = legendre(ell-1,cos_theta,'unnorm');
        previous_value = reshape( ...
            previous_values(absolute_order+1,:),1,[]);
    else
        previous_value = zeros(size(cos_theta));
    end
    associated_theta_derivative = ( ...
        ell*cos_theta.*associated_value- ...
        (ell+absolute_order)*previous_value)./sin_theta;
end

normalization = sqrt((2*ell+1)/(4*pi)*exp( ...
    gammaln(ell-absolute_order+1)-gammaln(ell+absolute_order+1)));
phase = exp(1i*absolute_order*phi);
Y_positive = normalization*associated_value.*phase;
theta_derivative_positive = ...
    normalization*associated_theta_derivative.*phase;

if emm<0
    negative_order_factor = (-1)^absolute_order;
    Y = negative_order_factor*conj(Y_positive);
    theta_derivative = negative_order_factor* ...
        conj(theta_derivative_positive);
else
    Y = Y_positive;
    theta_derivative = theta_derivative_positive;
end

phi_derivative = 1i*emm*Y;
theta_direction = [ ...
    cos_theta.*cos(phi); ...
    cos_theta.*sin(phi); ...
    -sin_theta];
phi_direction = [-sin(phi);cos(phi);zeros(size(phi))];
surface_gradient_Y = theta_direction.*theta_derivative+ ...
    phi_direction.*(phi_derivative./sin_theta);
end


function [nodes,weights] = gauss_legendre(number_of_nodes)

indices = 1:number_of_nodes-1;
off_diagonal = indices./sqrt(4*indices.^2-1);
jacobi_matrix = diag(off_diagonal,1)+diag(off_diagonal,-1);
[eigenvectors,eigenvalues] = eig(jacobi_matrix);
[nodes,permutation] = sort(diag(eigenvalues));
weights = 2*eigenvectors(1,permutation).^2;
nodes = nodes(:);
weights = weights(:);
end
