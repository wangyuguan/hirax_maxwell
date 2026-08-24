
clear
clc

test_directory = fileparts(mfilename('fullpath'));
fmm3dbie_directory = fullfile(test_directory,'..', ...
    'fmm3dbie-hirax-dev','matlab');
run(fullfile(fmm3dbie_directory,'startup.m'))

% Shared high-order unit-sphere discretization.
sphere_radius = 1;
number_of_cube_subdivisions = 2;
surface_order = 12;
patch_type = 1;
quadrature_tolerance = 1e-10;
gmres_tolerance = 1e-11;
required_solver_error = 1e-7;
required_evaluator_error = 2e-7;

S = geometries.sphere(sphere_radius,number_of_cube_subdivisions, ...
    [0;0;0],surface_order,patch_type);

fprintf('Maxwell unit-sphere analytic tests\n')
fprintf('  patches: %d, points: %d, order: %d\n', ...
    S.npatches,S.npts,surface_order)

%% PEC NRCCIE solver and near-corrected postprocessing


zk = 1.1;
alpha = 1;

pec_source = struct();
pec_source.r = [3.17;-0.03;3.15];
pec_source.edips = [1+0.30i;-0.40+0.20i;0.25-0.55i];
pec_source.hdips = [0.60-0.10i;-0.20+0.80i;0.40+0.30i];

[electric_incident,magnetic_incident] = em3d.incoming_sources( ...
    zk,pec_source,S,'electric and magnetic dipole');

pec_solver_options = struct();
pec_solver_options.rep = 'nrccie';
pec_solver_options.eps_gmres = gmres_tolerance;
pec_solver_options.maxit = 300;

[pec_densities,pec_gmres_history,pec_relative_residual] = ...
    em3d.pec.solver(S,electric_incident,magnetic_incident, ...
    quadrature_tolerance,zk,alpha,pec_solver_options);

far_targets = [ ...
     0.00,  0.17, -0.20,  0.12, -0.08,  0.13; ...
     0.00,  0.10,  0.15, -0.18, -0.12,  0.06; ...
     0.00, -0.11,  0.05,  0.06, -0.16, -0.14];

near_directions = [ ...
     1.00,  0.31, -0.44,  0.23, -0.62,  0.41; ...
     0.21, -0.81,  0.37,  0.91,  0.24, -0.55; ...
     0.17,  0.49,  0.82, -0.34,  0.71,  0.72];
near_directions = near_directions./vecnorm(near_directions,2,1);
near_radii = [0.990,0.985,0.980,0.975,0.970,0.965];
near_targets = near_directions.*near_radii;

pec_target_info = struct();
pec_target_info.r = [far_targets,near_targets];

pec_eval_options = struct('rep','nrccie');
[pec_electric,pec_magnetic] = em3d.pec.eval(S,pec_densities, ...
    pec_target_info,quadrature_tolerance,zk,alpha,pec_eval_options);
[pec_electric_incident,pec_magnetic_incident] = ...
    em3d.incoming_sources(zk,pec_source,pec_target_info, ...
    'electric and magnetic dipole');

nfar = size(far_targets,2);
far_indices = 1:nfar;
near_indices = nfar+(1:size(near_targets,2));

pec_far_error = relative_field_error( ...
    pec_electric(:,far_indices)+pec_electric_incident(:,far_indices), ...
    pec_magnetic(:,far_indices)+pec_magnetic_incident(:,far_indices), ...
    pec_electric_incident(:,far_indices), ...
    pec_magnetic_incident(:,far_indices));
pec_near_error = relative_field_error( ...
    pec_electric(:,near_indices)+pec_electric_incident(:,near_indices), ...
    pec_magnetic(:,near_indices)+pec_magnetic_incident(:,near_indices), ...
    pec_electric_incident(:,near_indices), ...
    pec_magnetic_incident(:,near_indices));

fprintf('\nPEC NRCCIE\n')
fprintf('  GMRES iterations: %d, relative residual: %.3e\n', ...
    numel(pec_gmres_history),pec_relative_residual)
fprintf('  far-target field error:  %.3e\n',pec_far_error)
fprintf('  near-target field error: %.3e\n',pec_near_error)

assert(isfinite(pec_relative_residual) && ...
    pec_relative_residual <= 10*gmres_tolerance, ...
    'NRCCIE: GMRES did not reach the requested tolerance.')
assert(isfinite(pec_far_error) && pec_far_error < required_solver_error, ...
    'NRCCIE: far-target field error %.3e is too large.',pec_far_error)
assert(isfinite(pec_near_error) && ...
    pec_near_error < required_evaluator_error, ...
    'NRCCIE: near-target field error %.3e is too large.',pec_near_error)

%% Isolated common-evaluator near-field regression
%
% For J=(1,0,0), rho=0 on the unit sphere,
%
%   S_k[1](x) = exp(i*k) j_0(k*r),  r < 1,
%   E = i*k*S_k[1]*(1,0,0),
%   H = (0, partial_z S_k[1], -partial_y S_k[1]).
%
% The final component specifically exercises
% H_z = partial_x S[J_y] - partial_y S[J_x].

evaluator_densities = complex(zeros(4,S.npts));
evaluator_densities(1,:) = 1;
evaluator_target_info = struct('r',near_targets);

[evaluator_electric,evaluator_magnetic] = em3d.pec.eval( ...
    S,evaluator_densities,evaluator_target_info,quadrature_tolerance, ...
    zk,alpha,pec_eval_options);

smooth_eval_options = struct('rep','nrccie','nonsmoothonly',true);
[smooth_electric,smooth_magnetic] = em3d.pec.eval( ...
    S,evaluator_densities,evaluator_target_info,quadrature_tolerance, ...
    zk,alpha,smooth_eval_options);

[evaluator_electric_exact,evaluator_magnetic_exact] = ...
    constant_current_sphere_field(near_targets,zk);

evaluator_electric_error = relative_array_error( ...
    evaluator_electric,evaluator_electric_exact);
evaluator_magnetic_error = relative_array_error( ...
    evaluator_magnetic,evaluator_magnetic_exact);
smooth_electric_error = relative_array_error( ...
    smooth_electric,evaluator_electric_exact);
smooth_magnetic_error = relative_array_error( ...
    smooth_magnetic,evaluator_magnetic_exact);

magnetic_scale = norm(evaluator_magnetic_exact,'fro');
magnetic_component_errors = vecnorm( ...
    evaluator_magnetic-evaluator_magnetic_exact,2,2)/magnetic_scale;

fprintf('\nCommon near-field evaluator\n')
fprintf('  corrected relative E error: %.3e\n',evaluator_electric_error)
fprintf('  corrected relative H error: %.3e\n',evaluator_magnetic_error)
fprintf('  corrected H component errors: [%.3e %.3e %.3e]\n', ...
    magnetic_component_errors)
fprintf('  smooth-only relative E/H errors: %.3e / %.3e\n', ...
    smooth_electric_error,smooth_magnetic_error)

assert(evaluator_electric_error < required_evaluator_error, ...
    'Common evaluator: near-corrected electric field is inaccurate.')
assert(evaluator_magnetic_error < required_evaluator_error, ...
    'Common evaluator: near-corrected magnetic field is inaccurate.')

%% Dielectric Muller transmission solver and postprocessing
%
% Manufacture an exact transmission solution from two point sources:
% a medium-0 source inside the sphere defines the radiating exterior field,
% and a medium-1 source outside the sphere defines the regular interior
% field. Their tangential jump supplies the incident boundary data.

omega = 1.1;
epsilon0 = 1.0;
mu0 = 1.0;
epsilon1 = 1.2;
mu1 = 1.0;
dielectric_parameters = complex([epsilon0,mu0,epsilon1,mu1]);

interior_source = struct();
interior_source.r = [0.11;0.00;0.37];
interior_source.edips = [1.00+0.20i;-0.35+0.15i;0.45-0.10i];
interior_source.hdips = [0.25-0.30i;0.70+0.10i;-0.20+0.40i];

exterior_source = struct();
exterior_source.r = [-3.50;3.10;5.10];
exterior_source.edips = interior_source.edips;
exterior_source.hdips = interior_source.hdips;

[exterior_trace_electric,exterior_trace_magnetic] = ...
    medium_dipole_fields(omega,epsilon0,mu0,interior_source,S);
[interior_trace_electric,interior_trace_magnetic] = ...
    medium_dipole_fields(omega,epsilon1,mu1,exterior_source,S);

transmission_electric_data = ...
    exterior_trace_electric-interior_trace_electric;
transmission_magnetic_data = ...
    exterior_trace_magnetic-interior_trace_magnetic;

dielectric_solver_options = struct();
dielectric_solver_options.rep = 'muller';
dielectric_solver_options.eps_gmres = gmres_tolerance;
dielectric_solver_options.maxit = 300;

[dielectric_densities,dielectric_gmres_history, ...
    dielectric_relative_residual] = em3d.dielectric.solver( ...
    S,transmission_electric_data,transmission_magnetic_data, ...
    quadrature_tolerance,omega,dielectric_parameters, ...
    dielectric_solver_options);

dielectric_interior_targets = [ ...
    -0.21,  0.15,  0.28, -0.33; ...
     0.18, -0.26,  0.11,  0.20; ...
    -0.14,  0.09, -0.24,  0.31];
dielectric_exterior_targets = [ ...
     2.10, -2.25,  1.80, -1.65; ...
     0.35,  0.80, -1.25, -0.55; ...
    -0.40,  0.30,  0.75, -1.20];
dielectric_target_info = struct();
dielectric_target_info.r = [ ...
    dielectric_interior_targets,dielectric_exterior_targets];

ninterior = size(dielectric_interior_targets,2);
interior_indices = 1:ninterior;
exterior_indices = ninterior+(1:size(dielectric_exterior_targets,2));
dielectric_eval_options = struct('rep','muller','in',interior_indices);

[dielectric_electric,dielectric_magnetic] = em3d.dielectric.eval( ...
    S,dielectric_densities,dielectric_target_info, ...
    quadrature_tolerance,omega,dielectric_parameters, ...
    dielectric_eval_options);

interior_exact_info = struct('r',dielectric_interior_targets);
exterior_exact_info = struct('r',dielectric_exterior_targets);
[interior_electric_exact,interior_magnetic_exact] = ...
    medium_dipole_fields(omega,epsilon1,mu1, ...
    exterior_source,interior_exact_info);
[exterior_electric_exact,exterior_magnetic_exact] = ...
    medium_dipole_fields(omega,epsilon0,mu0, ...
    interior_source,exterior_exact_info);

dielectric_interior_error = relative_field_error( ...
    dielectric_electric(:,interior_indices)-interior_electric_exact, ...
    dielectric_magnetic(:,interior_indices)-interior_magnetic_exact, ...
    interior_electric_exact,interior_magnetic_exact);
dielectric_exterior_error = relative_field_error( ...
    dielectric_electric(:,exterior_indices)-exterior_electric_exact, ...
    dielectric_magnetic(:,exterior_indices)-exterior_magnetic_exact, ...
    exterior_electric_exact,exterior_magnetic_exact);

fprintf('\nDielectric Muller\n')
fprintf('  GMRES iterations: %d, relative residual: %.3e\n', ...
    numel(dielectric_gmres_history),dielectric_relative_residual)
fprintf('  interior field error: %.3e\n',dielectric_interior_error)
fprintf('  exterior field error: %.3e\n',dielectric_exterior_error)

assert(isfinite(dielectric_relative_residual) && ...
    dielectric_relative_residual <= 20*gmres_tolerance, ...
    'Muller: GMRES did not reach the requested tolerance.')
assert(isfinite(dielectric_interior_error) && ...
    dielectric_interior_error < required_solver_error, ...
    'Muller: interior field error %.3e is too large.', ...
    dielectric_interior_error)
assert(isfinite(dielectric_exterior_error) && ...
    dielectric_exterior_error < required_solver_error, ...
    'Muller: exterior field error %.3e is too large.', ...
    dielectric_exterior_error)

fprintf('\nPASS: every implemented MATLAB Maxwell solver and the corrected ')
fprintf('near-field evaluator agree with analytic sphere solutions.\n')


function error_value = relative_array_error(computed,exact)
error_value = norm(computed-exact,'fro')/norm(exact,'fro');
end


function error_value = relative_field_error( ...
    electric_difference,magnetic_difference, ...
    electric_reference,magnetic_reference)
error_value = norm([electric_difference;magnetic_difference],'fro')/ ...
    norm([electric_reference;magnetic_reference],'fro');
end


function [electric,magnetic] = constant_current_sphere_field(targets,zk)
r = vecnorm(targets,2,1);
z = zk*r;
spherical_bessel_j0 = sin(z)./z;
spherical_bessel_j0_derivative = (z.*cos(z)-sin(z))./z.^2;
scalar_single_layer = exp(1i*zk)*spherical_bessel_j0;
radial_derivative = exp(1i*zk)*zk*spherical_bessel_j0_derivative;
gradient_single_layer = targets./r.*radial_derivative;

electric = complex(zeros(3,size(targets,2)));
electric(1,:) = 1i*zk*scalar_single_layer;
magnetic = complex(zeros(3,size(targets,2)));
magnetic(2,:) = gradient_single_layer(3,:);
magnetic(3,:) = -gradient_single_layer(2,:);
end


function [electric,magnetic] = medium_dipole_fields( ...
    omega,epsilon,mu,source,target_info)
% Match fieldsEDomega/fieldsMDomega using the normalized MATLAB dipoles.
wave_number = omega*sqrt(epsilon*mu);
number_of_targets = size(target_info.r,2);
electric = complex(zeros(3,number_of_targets));
magnetic = complex(zeros(3,number_of_targets));

if isfield(source,'edips')
    electric_source = struct('r',source.r,'edips',source.edips);
    [electric_ed,magnetic_ed] = em3d.incoming_sources( ...
        wave_number,electric_source,target_info,'electric dipole');
    electric = electric-1i/(omega*epsilon)*electric_ed;
    magnetic = magnetic-1i/wave_number*magnetic_ed;
end

if isfield(source,'hdips')
    magnetic_source = struct('r',source.r,'hdips',source.hdips);
    [electric_hd,magnetic_hd] = em3d.incoming_sources( ...
        wave_number,magnetic_source,target_info,'magnetic dipole');
    electric = electric-1i/wave_number*electric_hd;
    magnetic = magnetic-1i/(omega*mu)*magnetic_hd;
end
end
