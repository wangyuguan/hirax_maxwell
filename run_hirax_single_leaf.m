
clear
clc

run('../fmm3dbie-hirax-dev/matlab/startup.m')
run('../chunkie/startup.m')
addpath('src')

dish_radius = 69.25;
thickness = 0.01*dish_radius;
surface_orders = [6 8 10 12 14 16];

chunkie_order = 20;
chunkie_n0 = 3;
chunkie_nchs = 3;
chunkie_newton_iterations = 30;
rim_width = 0.028128271246*dish_radius;
cap_collar_width = rim_width;
cap_mesh_spacing_center = 0.2*dish_radius;
cap_mesh_spacing_side = cap_mesh_spacing_center;
cap_mesh_side_start = 0.62;
outline_refinement = 1;
wall_profile_refinement = 3;

wavelength = 10*dish_radius;
zk = 2*pi/wavelength;
source_gap = dish_radius/30;
source_height = 0;
outline_samples_per_panel = 2001;

alpha = 1.0;
eps_quad = 1e-7;
eps_gmres = 1e-8;
maxit = 300;

p0 = dish_radius^3*[1;1i;0];
test_targets = dish_radius*[ ...
     0.0,  0.0,  2.0,  0.0, -1.6,  1.2; ...
     0.0,  0.0,  0.0,  2.0,  0.7, -1.3; ...
     2.0, -2.0,  0.4, -0.3,  0.8, -1.1];

data_file = 'run_hirax_single_leaf_data.mat';
log_file = 'run_hirax_single_leaf.log';
level_directory = 'data';

if ~isfolder(level_directory)
    mkdir(level_directory)
end

diary(log_file)
diary on
diary_cleanup = onCleanup(@()diary('off'));

% Construct the fixed outline once and locate its actual rightmost point.
% This uses the Chunkie-smoothed blue curve, not the original analytic
% circle-ellipse join. The same master outline is used at every order.
reference_options = struct();
reference_options.norder = surface_orders(1);
reference_options.chunkie_order = chunkie_order;
reference_options.chunkie_n0 = chunkie_n0;
reference_options.chunkie_nchs = chunkie_nchs;
reference_options.chunkie_newton_iterations = ...
    chunkie_newton_iterations;
reference_options.rim_width = rim_width;
reference_options.cap_collar_width = cap_collar_width;
reference_options.cap_mesh_spacing_center = cap_mesh_spacing_center;
reference_options.cap_mesh_spacing_side = cap_mesh_spacing_side;
reference_options.cap_mesh_side_start = cap_mesh_side_start;
reference_options.outline_refinement = outline_refinement;
reference_options.wall_profile_refinement = wall_profile_refinement;

reference_timer = tic;
[reference_surface,reference_parts] = ...
    hirax_chunkie_dish_plate_surfer(thickness,reference_options);
reference_geometry_time = toc(reference_timer);
source_anchor = rightmost_outline_point( ...
    reference_parts,outline_samples_per_panel);
x0 = source_anchor+[source_gap;0;source_height];

source_info = struct();
source_info.r = x0;
source_info.edips = -p0;

target_info = struct();
target_info.r = test_targets;

[einc_probe,hinc_probe] = em3d.incoming_sources( ...
    zk,source_info,target_info,'electric dipole');

settings = struct();
settings.thickness = thickness;
settings.surface_orders = surface_orders;
settings.chunkie_order = chunkie_order;
settings.chunkie_n0 = chunkie_n0;
settings.chunkie_nchs = chunkie_nchs;
settings.chunkie_newton_iterations = chunkie_newton_iterations;
settings.rim_width = rim_width;
settings.cap_collar_width = cap_collar_width;
settings.cap_mesh_spacing_center = cap_mesh_spacing_center;
settings.cap_mesh_spacing_side = cap_mesh_spacing_side;
settings.cap_mesh_side_start = cap_mesh_side_start;
settings.outline_refinement = outline_refinement;
settings.wall_profile_refinement = wall_profile_refinement;
settings.dish_radius = dish_radius;
settings.wavelength = wavelength;
settings.zk = zk;
settings.source_gap = source_gap;
settings.source_height = source_height;
settings.source_anchor = source_anchor;
settings.x0 = x0;
settings.p0 = p0;
settings.alpha = alpha;
settings.eps_quad = eps_quad;
settings.eps_gmres = eps_gmres;
settings.maxit = maxit;
settings.test_targets = test_targets;
settings.einc_probe = einc_probe;
settings.hinc_probe = hinc_probe;

number_of_orders = numel(surface_orders);
levels = cell(number_of_orders,1);

fprintf('HIRAX Chunkie dish close-source NRCCIE raw-data run\n')
fprintf('surface orders %s, thickness %.8g\n', ...
    mat2str(surface_orders),thickness)
fprintf(['dish radius %.8g, wavelength %.8g, R/lambda %.8g, ' ...
    'zk %.8g\n'],dish_radius,wavelength,dish_radius/wavelength,zk)
fprintf(['source anchor (%.12g, %.12g, %.12g), gap R/30 = %.12g\n' ...
    'source point  (%.12g, %.12g, %.12g)\n'], ...
    source_anchor,source_gap,x0)
fprintf(['reference geometry %.1f s, eps_quad %.1e, ' ...
    'eps_gmres %.1e\n'],reference_geometry_time,eps_quad,eps_gmres)
fprintf(['adaptive cap spacing: center %.6g, side %.6g, ' ...
    'transition starts at |x|/xmax = %.3g\n'], ...
    cap_mesh_spacing_center,cap_mesh_spacing_side,cap_mesh_side_start)
fprintf(['uniform outline subdivision: %d, ' ...
    'one cap collar layer, wall profile panels per half: %d\n'], ...
    outline_refinement,wall_profile_refinement)
fprintf('level directory: %s\n\n',level_directory)

for order_id = 1:number_of_orders
    surface_order = surface_orders(order_id);
    level_file = sprintf('%s/run_hirax_single_leaf_order%d.mat', ...
        level_directory,surface_order);

    case_settings = rmfield(settings,'surface_orders');
    case_settings.surface_order = surface_order;

    fprintf('\n========== surface order %d ==========%s', ...
        surface_order,newline)

    try
        geometry_options = struct();
        geometry_options.norder = surface_order;
        geometry_options.chunkie_order = chunkie_order;
        geometry_options.chunkie_n0 = chunkie_n0;
        geometry_options.chunkie_nchs = chunkie_nchs;
        geometry_options.chunkie_newton_iterations = ...
            chunkie_newton_iterations;
        geometry_options.rim_width = rim_width;
        geometry_options.cap_collar_width = cap_collar_width;
        geometry_options.cap_mesh_spacing_center = ...
            cap_mesh_spacing_center;
        geometry_options.cap_mesh_spacing_side = cap_mesh_spacing_side;
        geometry_options.cap_mesh_side_start = cap_mesh_side_start;
        geometry_options.outline_refinement = outline_refinement;
        geometry_options.wall_profile_refinement = ...
            wall_profile_refinement;

        if surface_order==surface_orders(1) && ...
                exist('reference_surface','var')
            S = reference_surface;
            plate_parts = reference_parts;
            geometry_time = reference_geometry_time;
            clear reference_surface reference_parts
        else
            geometry_timer = tic;
            [S,plate_parts] = hirax_chunkie_dish_plate_surfer( ...
                thickness,geometry_options);
            geometry_time = toc(geometry_timer);
        end

        surface_weights = S.wts(:).';
        surface_area = sum(surface_weights);
        signed_volume = sum(sum(S.r.*S.n,1).*surface_weights)/3;
        sampled_source_distance = min(vecnorm(S.r-x0,2,1));
        probe_distances = zeros(1,size(test_targets,2));
        for target_id = 1:size(test_targets,2)
            probe_distances(target_id) = min(vecnorm( ...
                S.r-test_targets(:,target_id),2,1));
        end

        fprintf(['requested source gap %.12g, nearest surface node ' ...
            '%.12g\n'],source_gap,sampled_source_distance)
        fprintf(['mesh: %d core triangles per face, %d outline panels, ' ...
            'one cap collar layer, %d wall profile panels per half\n'], ...
            plate_parts.number_of_core_triangles, ...
            plate_parts.number_of_outline_panels, ...
            plate_parts.number_of_wall_profile_panels)

        [einc,hinc] = em3d.incoming_sources( ...
            zk,source_info,S,'electric dipole');

        solver_options = struct();
        solver_options.rep = 'nrccie';
        solver_options.eps_gmres = eps_gmres;
        solver_options.maxit = maxit;

        solve_timer = tic;
        [densities,gmres_history,relative_residual,Q,solver_timing] = ...
            em3d.pec.solver( ...
            S,einc,hinc,eps_quad,zk,alpha,solver_options);
        measured_solve_time = toc(solve_timer);
        iterations = numel(gmres_history);
        converged = isfinite(relative_residual) && ...
            relative_residual<=eps_gmres;

        if isfield(Q,'nnz')
            quadrature_nnz = Q.nnz;
        else
            quadrature_nnz = length(Q.col_ind);
        end
        quadrature_nquad = Q.nquad;

        escat = complex(nan(3,size(test_targets,2)));
        hscat = complex(nan(3,size(test_targets,2)));
        evaluation_time = NaN;
        evaluation_finished = false;
        evaluation_error_identifier = '';
        evaluation_error_message = '';

        if converged
            evaluator_options = struct();
            evaluator_options.rep = 'nrccie';
            evaluator_options.nonsmoothonly = true;
            try
                evaluation_timer = tic;
                [escat,hscat] = em3d.pec.eval( ...
                    S,densities,target_info,eps_quad,zk,alpha, ...
                    evaluator_options);
                evaluation_time = toc(evaluation_timer);
                evaluation_finished = all(isfinite([escat(:);hscat(:)]));
            catch evaluation_error
                evaluation_error_identifier = evaluation_error.identifier;
                evaluation_error_message = evaluation_error.message;
            end
        end

        if ~converged
            status = 'gmres_not_converged';
        elseif ~evaluation_finished
            status = 'evaluation_failed';
        else
            status = 'complete';
        end

        level = struct();
        level.finished = true;
        level.status = status;
        level.surface_order = surface_order;
        level.npatches = S.npatches;
        level.npts = S.npts;
        level.iptype = S.iptype;
        level.number_of_core_triangles_per_face = ...
            plate_parts.number_of_core_triangles;
        level.number_of_outline_panels = ...
            plate_parts.number_of_outline_panels;
        level.number_of_cap_collar_layers = ...
            plate_parts.number_of_cap_collar_layers;
        level.number_of_wall_profile_panels_per_half = ...
            plate_parts.number_of_wall_profile_panels;
        level.surface_area = surface_area;
        level.signed_volume = signed_volume;
        level.max_position_mismatch = plate_parts.max_position_mismatch;
        level.max_normal_mismatch_degrees = ...
            plate_parts.max_normal_mismatch_degrees;
        level.source_anchor = source_anchor;
        level.source_point = x0;
        level.requested_source_gap = source_gap;
        level.source_distance = sampled_source_distance;
        level.sampled_source_distance = sampled_source_distance;
        level.probe_distances = probe_distances;
        level.iterations = iterations;
        level.gmres_history = gmres_history;
        level.relative_residual = relative_residual;
        level.converged = converged;
        level.density_norm = norm(densities,'fro');
        level.quadrature_nnz = quadrature_nnz;
        level.quadrature_nquad = quadrature_nquad;
        level.geometry_time = geometry_time;
        level.quadrature_time = solver_timing.quadrature_time;
        level.oversampling_time = solver_timing.oversampling_time;
        level.gmres_time = solver_timing.gmres_time;
        level.solver_total_time = solver_timing.total_time;
        level.measured_solve_time = measured_solve_time;
        level.evaluation_time = evaluation_time;
        level.evaluation_finished = evaluation_finished;
        level.evaluation_error_identifier = evaluation_error_identifier;
        level.evaluation_error_message = evaluation_error_message;
        level.escat = escat;
        level.hscat = hscat;
        level.etotal = einc_probe+escat;
        level.htotal = hinc_probe+hscat;
        levels{order_id} = level;

        save(level_file,'case_settings','level','S','densities', ...
            'einc','hinc','-v7.3')
        save(data_file,'settings','levels','-v7.3')

        fprintf(['order %d saved: status %s, %d patches, %d nodes\n' ...
            'GMRES %d iterations, residual %.3e, solve %.1f s\n' ...
            'probe evaluation %.1f s, file %s\n'], ...
            surface_order,status,S.npatches,S.npts,iterations, ...
            relative_residual,measured_solve_time,evaluation_time,level_file)

        clear S densities Q einc hinc
    catch run_error
        level = struct();
        level.finished = false;
        level.status = 'run_error';
        level.surface_order = surface_order;
        level.error_identifier = run_error.identifier;
        level.error_message = run_error.message;
        level.error_report = getReport(run_error,'extended', ...
            'hyperlinks','off');
        levels{order_id} = level;
        save(level_file,'case_settings','level','-v7.3')
        save(data_file,'settings','levels','-v7.3')
        fprintf(2,'order %d failed; diagnostic saved in %s\n%s\n', ...
            surface_order,level_file,level.error_report)
        rethrow(run_error)
    end
end

save(data_file,'settings','levels','-v7.3')
fprintf('\nClose-source raw-data run finished. Saved %s\n',data_file)


function point = rightmost_outline_point(parts,points_per_panel)
q = linspace(0,1,points_per_panel);
polynomials = lege.pols(2*q-1,parts.chunkie_order-1);
point = [-inf;NaN;0];
for panel = 1:parts.number_of_outline_panels
    coefficients = parts.outline_position_coefficients(:,:,panel);
    values = coefficients*reshape( ...
        polynomials,parts.chunkie_order,[]);
    [candidate_x,index] = max(values(1,:));
    if candidate_x>point(1)
        point = [values(:,index);0];
    end
end
end
