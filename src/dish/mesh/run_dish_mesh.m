function result = run_dish_mesh(stem, varargin)
%RUN_DISH_MESH High-order mesh of a HIRAX parabolic dish, graded at the rim.
%
%   result = run_dish_mesh('dish_fd021_m')
%   result = run_dish_mesh('dish_fd021_t1', 'order', 8, 'h_max', 5, ...
%                          'h_rim', 0.5, 'target', 3.4)
%
% Meshes ../cad/<stem>.step with step_mesher_cpp_matlab and exports the
% fmm3dbie formats. The CAD is built by ../cad/make_dish_geo.py, which also
% writes ../cad/<stem>.json holding the exact analytic geometry; the areas in
% there are what this driver validates against, rather than any CAD-reported
% number.
%
% Dependencies already on the MATLAB path are reused. Otherwise repositories
% are sought under NRCCIE_DEPENDENCY_ROOT (the parent of this project by
% default): fmm3dbie-hirax-dev and step_mesher_cpp_matlab. STEPMESHER_ROOT can
% explicitly select a different step_mesher_cpp_matlab checkout.
%
% The size field is a distance field measured from the rim face, the only face
% small enough to be picked up by refine_face_area_max: elements are h_rim on
% and next to the rim and coarsen to the far-field size away from it.
%
% 'target' is the far-field size handed to gmsh. gmsh treats a mesh size as a
% target, not a bound, and the realised longest edge overshoots it, so 'target'
% is set below 'h_max' and the driver then checks the realised element
% diameters against 'h_max' and reports whether the bound actually holds. Both
% the flat scaffold triangles and the curved high-order patches are measured:
% the patch diameter is the largest distance between any two of its projected
% nodes, which is the honest element diameter of the surface that gets solved on.

p = inputParser;
p.addParameter('order', 8, @(x) isscalar(x) && x > 0);
p.addParameter('h_max', 0.05, @isscalar);       % required bound on element diameter
p.addParameter('h_rim', 0.01, @isscalar);       % element size at the rim
p.addParameter('target', 0.05, @(x) isempty(x) || isscalar(x));      % far-field size handed to gmsh
p.addParameter('dist_min', 0.02, @isscalar);    % hold h_rim within this of the rim
p.addParameter('dist_max', 0.25, @isscalar);   % reach the far-field size here
p.addParameter('profile', 'rigid', @ischar);
p.addParameter('edge_gll_order', 16, @isscalar);
p.addParameter('mesh_fraction', [], @(x) isempty(x) || isscalar(x));
p.addParameter('tag', '', @ischar);        % suffix for the output dir/filenames
p.addParameter('plots', true, @islogical);
p.addParameter('export', true, @islogical);
p.parse(varargin{:});
cfg = p.Results;

%% Paths
thisDir  = fileparts(mfilename('fullpath'));
dishRoot = fileparts(thisDir);
stepFile = fullfile(dishRoot, 'cad', [stem '.step']);
jsonFile = fullfile(dishRoot, 'cad', [stem '.json']);
assert(isfile(stepFile), 'missing STEP file %s', stepFile);
assert(isfile(jsonFile), 'missing CAD record %s', jsonFile);

projectRoot = fileparts(fileparts(dishRoot));
dependencyRoot = getenv('NRCCIE_DEPENDENCY_ROOT');
if isempty(dependencyRoot), dependencyRoot = fileparts(projectRoot); end

mesherRoot = getenv('STEPMESHER_ROOT');
if ~isempty(mesherRoot) || isempty(which('stepmesher.mesh_step'))
    if isempty(mesherRoot)
        mesherRoot = fullfile(dependencyRoot, 'step_mesher_cpp_matlab');
    end
    mesherPath = fullfile(mesherRoot, 'matlab');
    assert(isfolder(mesherPath), ...
        'missing mesher MATLAB folder %s; set STEPMESHER_ROOT', mesherPath);
    addpath(mesherPath);
end
assert(~isempty(which('stepmesher.mesh_step')), ...
    'stepmesher.mesh_step is unavailable; check STEPMESHER_ROOT');

% fmm3dbie supplies the koorn package the RV nodes and exports need.
if isempty(which('koorn.rv_nodes')) || isempty(which('surfer'))
    fmm3dbiePath = fullfile(dependencyRoot, 'fmm3dbie-hirax-dev', 'matlab');
    startupFile = fullfile(fmm3dbiePath, 'startup.m');
    assert(isfile(startupFile), ...
        'missing fmm3dbie startup %s; set NRCCIE_DEPENDENCY_ROOT', startupFile);
    run(startupFile);
end

cad = jsondecode(fileread(jsonFile));

outName = [stem cfg.tag];
outDir  = fullfile(thisDir, 'out', outName);
if ~isfolder(outDir), mkdir(outDir); end

%% Options
opts = struct();
opts.order            = cfg.order;
opts.refinement_level = 0;
opts.edge_gll_order   = cfg.edge_gll_order;
opts.mesh_profile     = cfg.profile;
% h_rim <= 0 leaves the refinement field out entirely: step_mesher then keeps
% the uniform mesh_fraction*bbox_diag size everywhere, which is the
% "no size restriction" case.
if cfg.h_rim > 0
    opts.refine_size     = cfg.h_rim;
    % refine_size_max = 0 makes step_mesher fall back to the uniform
    % mesh_fraction*bbox_diag size as the far field, i.e. grade the rim but
    % leave the bulk unconstrained.
    if isempty(cfg.target) || cfg.target <= 0
        opts.refine_size_max = 0;
    else
        opts.refine_size_max = cfg.target;
    end
    opts.refine_dist_min = cfg.dist_min;
    opts.refine_dist_max = cfg.dist_max;
else
    opts.refine_size = 0;
end
% The rim is the only face small enough to be caught: it is about pi*D*t
% against ~D^2 for each paraboloid, so any cutoff between those two selects
% exactly the rim without needing its tag. Taking their geometric mean keeps
% that true whatever length unit the CAD carries.
opts.refine_face_area_max = sqrt(cad.area_rim_exact * ...
    min(cad.area_inner_exact, cad.area_outer_exact));
% Irrelevant while the refinement field is active (it supplies both the min and
% max size), but set so the uniform fallback would also land on the target.
% With the refinement field active this is irrelevant (the field supplies both
% the min and the max size); with it off this is the only thing setting the
% element size, so an explicit mesh_fraction overrides the target-derived one.
bboxDiag = hypot(hypot(2*cad.rim_outer_r, 2*cad.rim_outer_r), ...
                 cad.rim_inner_z - cad.vertex_z_outer);
if ~isempty(cfg.mesh_fraction)
    opts.mesh_fraction = cfg.mesh_fraction;
elseif ~isempty(cfg.target) && cfg.target > 0
    opts.mesh_fraction = cfg.target / bboxDiag;
else
    opts.mesh_fraction = 0.10;    % step_mesher's own default
end

fprintf('\n=== %s ===\n', stem);
fprintf('f/D %.4g   D %.4g %s   f %.4g %s   wall %.4g %s\n', ...
    cad.f_ratio, cad.diameter, cad.units, cad.focal_length, cad.units, ...
    cad.wall_thickness, cad.units);
if isempty(cfg.target) || cfg.target <= 0
    targetStr = 'uniform (mesh_fraction)';
else
    targetStr = sprintf('%.4g', cfg.target);
end
if isinf(cfg.h_max)
    boundStr = 'none';
else
    boundStr = sprintf('%.4g', cfg.h_max);
end
fprintf('order %d   h_rim %.4g   far-field %s   bound h_max %s   (%s)\n', ...
    cfg.order, cfg.h_rim, targetStr, boundStr, cad.units);
if cfg.h_rim > 0
    farField = opts.refine_size_max;
    if farField <= 0, farField = opts.mesh_fraction * bboxDiag; end
    fprintf('grading: %.4g within %.4g of the rim, %.4g beyond %.4g (%s)\n\n', ...
        cfg.h_rim, cfg.dist_min, farField, cfg.dist_max, cad.units);
else
    fprintf('grading: none; uniform mesh_fraction %.6g -> %.4g %s\n\n', ...
        opts.mesh_fraction, opts.mesh_fraction*bboxDiag, cad.units);
end

%% Mesh
tm = tic;
mesh = stepmesher.mesh_step(stepFile, opts);
tmesh = toc(tm);
fprintf('meshed in %.1f s\n', tmesh);
fprintf('  patches           %d\n', mesh.npatches);
fprintf('  nodes per patch   %d\n', mesh.nodes_per_patch);
fprintf('  total nodes       %d\n', size(mesh.xyz, 2));
fprintf('  refined faces     %d\n', mesh.meta.refine_face_count);
fprintf('  continuity breaks %d\n', mesh.meta.continuity_break_count);
fprintf('  fragmented faces  %d\n', mesh.meta.fragmented_face_count);
fprintf('  projection fails  %d\n', mesh.meta.projection_failures);
fprintf('  max proj distance %.6e %s\n', mesh.meta.max_projection_distance, cad.units);
fprintf('  edge fallbacks    %d\n\n', mesh.meta.edge_fallback_points);

%% Element diameters
tri = mesh.scaffold_triangles;
X   = mesh.scaffold_nodes;
P1 = X(:, tri(1, :));  P2 = X(:, tri(2, :));  P3 = X(:, tri(3, :));
h_flat = max([vecnorm(P2 - P1); vecnorm(P3 - P2); vecnorm(P1 - P3)], [], 1);

% True diameter of each curved patch: largest distance between any two of its
% projected nodes. For a triangle the diameter is its longest edge, but a
% curved high-order patch can bulge, so measure all pairs rather than assume.
npp = mesh.nodes_per_patch;
R3  = reshape(mesh.xyz, 3, npp, mesh.npatches);
h_patch = zeros(1, mesh.npatches);
for a = 1:npp-1
    d = vecnorm(R3(:, a+1:npp, :) - R3(:, a, :), 2, 1);
    h_patch = max(h_patch, max(reshape(d, npp-a, mesh.npatches), [], 1));
end

% Aspect ratio per scaffold triangle: longest edge / shortest edge (1 for
% equilateral). A sliver spanning the 0.01 thick rim band with a long
% circumferential edge shows up here directly.
eL = [vecnorm(P2 - P1); vecnorm(P3 - P2); vecnorm(P1 - P3)];
aspect = max(eL, [], 1) ./ min(eL, [], 1);

%% Which face is the rim
faceTags = unique(mesh.scaffold_face_tags);
triArea  = 0.5 * vecnorm(cross(P2 - P1, P3 - P1), 2, 1);
faceArea = arrayfun(@(t) sum(triArea(mesh.scaffold_face_tags == t)), faceTags);
[~, irim] = min(abs(faceArea - cad.area_rim_exact));
rimTag = faceTags(irim);
onRim  = mesh.scaffold_face_tags == rimTag;

% Distance of each scaffold triangle from the rim circle, for the grading report.
C   = (P1 + P2 + P3) / 3;
rho = hypot(C(1, :), C(2, :));
dRim = hypot(rho - cad.rim_inner_r, C(3, :) - cad.rim_inner_z);

fprintf('element diameter (%s)\n', cad.units);
fprintf('  flat scaffold triangles   max %.4f   median %.4f   min %.4f\n', ...
    max(h_flat), median(h_flat), min(h_flat));
fprintf('  curved order-%d patches    max %.4f   median %.4f   min %.4f\n', ...
    cfg.order, max(h_patch), median(h_patch), min(h_patch));
fprintf('  bound h_max = %.4g %s:  %d / %d patches exceed it\n', cfg.h_max, cad.units, ...
    sum(h_patch > cfg.h_max), mesh.npatches);
if isinf(cfg.h_max)
    fprintf('  no bound requested (h_max = inf)\n\n');
elseif max(h_patch) <= cfg.h_max
    fprintf('  BOUND HOLDS (max patch diameter %.5g %s)\n\n', max(h_patch), cad.units);
else
    fprintf('  *** max patch diameter %.5g %s EXCEEDS h_max by %.5g ***\n\n', ...
        max(h_patch), cad.units, max(h_patch) - cfg.h_max);
end

fprintf('aspect ratio (longest edge / shortest edge, 1 = equilateral)\n');
fprintf('  whole mesh    median %.2f   p99 %.2f   max %.2f\n', ...
    median(aspect), prctile(aspect, 99), max(aspect));
fprintf('  rim face      median %.2f   p99 %.2f   max %.2f\n\n', ...
    median(aspect(onRim)), prctile(aspect(onRim), 99), max(aspect(onRim)));

fprintf('grading, by distance from the rim circle\n');
fprintf('  on the rim face (tag %d)  n = %6d   median %.4f   max %.4f\n', ...
    rimTag, sum(onRim), median(h_flat(onRim)), max(h_flat(onRim)));
edges = cfg.dist_max * [0 0.04 0.08 0.2 0.4 1 2 4 8 inf];
for k = 1:numel(edges)-1
    m = dRim >= edges(k) & dRim < edges(k+1) & ~onRim;
    if any(m)
        fprintf('  d in [%5.4g,%5.4g)        n = %6d   median %.4f   max %.4f\n', ...
            edges(k), edges(k+1), sum(m), median(h_flat(m)), max(h_flat(m)));
    end
end
fprintf('\n');

%% fmm3dbie conversion and validation
[srcvals, norders, iptype, ~] = stepmesher.to_srcvals(mesh); %#ok<ASGLU>
S = stepmesher.to_surfer(mesh);
stats = stepmesher.validate_area_flux(mesh);

relExact = abs(stats.area - cad.area_total_exact) / cad.area_total_exact;
fprintf('validation\n');
fprintf('  surfer            npatches %d   npts %d\n', S.npatches, S.npts);
fprintf('  quadrature area   %.10e %s^2\n', stats.area, cad.units);
fprintf('  exact CAD area    %.10e %s^2   (analytic, from %s.json)\n', ...
    cad.area_total_exact, cad.units, stem);
fprintf('  rel area error    %.6e   (vs exact)\n', relExact);
fprintf('  rel area error    %.6e   (vs cad_area reported by OCC)\n', ...
    stats.relative_area_error);
fprintf('  normal flux       [%.4e %.4e %.4e]  (0 for a closed surface)\n\n', ...
    stats.normal_flux(1), stats.normal_flux(2), stats.normal_flux(3));

%% Export
result = struct('stem', stem, 'out_name', outName, 'cad', cad, 'opts', opts, 'mesh_seconds', tmesh, ...
    'npatches', mesh.npatches, 'total_nodes', size(mesh.xyz, 2), ...
    'h_flat_max', max(h_flat), 'h_patch_max', max(h_patch), ...
    'h_patch_median', median(h_patch), 'h_bound', cfg.h_max, ...
    'aspect_median', median(aspect), 'aspect_p99', prctile(aspect, 99), ...
    'aspect_max', max(aspect), 'aspect_rim_median', median(aspect(onRim)), ...
    'aspect_rim_max', max(aspect(onRim)), ...
    'bound_holds', max(h_patch) <= cfg.h_max, ...
    'area', stats.area, 'rel_area_error_exact', relExact, ...
    'normal_flux', stats.normal_flux, 'rim_face_tag', rimTag);

if cfg.export
    xyzwFile  = fullfile(outDir, sprintf('%s_order%d_xyzw.txt', outName, cfg.order));
    go3File   = fullfile(outDir, sprintf('%s_order%d.go3', outName, cfg.order));
    gidFile   = fullfile(outDir, sprintf('%s_scaffold.gidmsh', outName));
    pxyzw = stepmesher.export_xyzw(mesh, xyzwFile);
    pgo3  = stepmesher.export_go3(mesh, go3File);
    pgid  = stepmesher.export_gidmsh(mesh, gidFile);
    fprintf('export\n');
    fprintf('  %s  (%d nodes)\n', xyzwFile, pxyzw.total_nodes);
    fprintf('  %s  (%d nodes)\n', go3File, pgo3.total_nodes);
    fprintf('  %s  (%d triangles)\n\n', gidFile, pgid.num_triangles);
    result.xyzw_file = xyzwFile;
    result.go3_file = go3File;
    result.gidmsh_file = gidFile;
    save(fullfile(outDir, sprintf('%s_order%d_report.mat', outName, cfg.order)), ...
        '-struct', 'result');
end

%% Figures
if cfg.plots
    figDir = fullfile(outDir, 'figures');
    if ~isfolder(figDir), mkdir(figDir); end

    fg = figure('Visible', 'off', 'Position', [100 100 900 700]);
    trisurf(tri.', X(1, :).', X(2, :).', X(3, :).', h_flat.', ...
        'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.1);
    axis equal; colormap(parula); cb = colorbar; cb.Label.String = sprintf('element diameter (%s)', cad.units);
    view(-35, 22); title(sprintf('%s scaffold, graded at the rim', ...
        strrep(stem, '_', '\_')));
    xlabel(sprintf('x (%s)', cad.units)); ylabel(sprintf('y (%s)', cad.units));
    zlabel(sprintf('z (%s)', cad.units));
    exportgraphics(fg, fullfile(figDir, [outName '_scaffold.png']), 'Resolution', 150);
    close(fg);

    % Rim close-up: the graded band is only visible zoomed in.
    fg = figure('Visible', 'off', 'Position', [100 100 900 700]);
    zoom = 0.12 * cad.aperture_radius;
    keep = rho > cad.aperture_radius - zoom;
    trisurf(tri(:, keep).', X(1, :).', X(2, :).', X(3, :).', h_flat(keep).', ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.25);
    axis equal; cb = colorbar; cb.Label.String = sprintf('element diameter (%s)', cad.units);
    view(-50, 12);
    xlim([cad.aperture_radius - 1.15*zoom, cad.aperture_radius + 0.15*zoom]);
    ylim(1.15*zoom*[-1 1]);
    title(sprintf('%s rim grading (outer %.3g %s)', strrep(stem, '_', '\_'), ...
        zoom, cad.units));
    xlabel(sprintf('x (%s)', cad.units)); ylabel(sprintf('y (%s)', cad.units));
    zlabel(sprintf('z (%s)', cad.units));
    exportgraphics(fg, fullfile(figDir, [outName '_rim.png']), 'Resolution', 150);
    close(fg);

    fg = figure('Visible', 'off', 'Position', [100 100 800 550]);
    semilogy(dRim, h_flat, '.', 'MarkerSize', 2); hold on
    if isfinite(cfg.h_max)
        yline(cfg.h_max, 'r-', 'h\_max', 'LineWidth', 1.2);
    end
    if cfg.h_rim > 0
        yline(cfg.h_rim, 'b--', 'h\_rim', 'LineWidth', 1.2);
    end
    xline(cfg.dist_min, 'k:'); xline(cfg.dist_max, 'k:');
    xlabel(sprintf('distance from the rim circle (%s)', cad.units));
    ylabel(sprintf('element diameter (%s)', cad.units));
    title(sprintf('%s grading profile', strrep(stem, '_', '\_')));
    xlim([0 4*cfg.dist_max]); grid on
    exportgraphics(fg, fullfile(figDir, [outName '_grading.png']), 'Resolution', 150);
    close(fg);

    fprintf('figures in %s\n\n', figDir);
end
end
