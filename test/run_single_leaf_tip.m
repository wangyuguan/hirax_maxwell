% Self-convergence of one distinct-tip leaf excited from a neighboring leaf.
clear
clc

test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(fullfile(root,'src'))
clear pth dir

%% Run settings
settings.geometry_revision = leaf_geometry_revision();
settings.geometry_options = leaf_geometry_options();
settings.leaf_id = 1;                       % TL/TR/BL/BR ordering
settings.surface_orders = [4 6 8 10];
settings.integration_orders = [16 20];
settings.thickness = .6925;                 % mm
settings.leaf_radius = 69.25;               % mm
settings.zk = 2*pi/(10*settings.leaf_radius);
settings.alpha = 1;
settings.eps_quad = 1e-11;
settings.eps_fmm = 1e-9;
settings.eps_gmres = 1e-8;
settings.gmres_restart = 100;
settings.gmres_maxit = 10;
settings.quad_batch_size = 2000;
settings.outline_samples_per_panel = 17;
settings.electric_dipole = -settings.leaf_radius^3*[1;1i;0];
settings.geometry_only = strcmp( ...
    getenv('HIRAX_SINGLE_LEAF_GEOMETRY_ONLY'),'1');

order_override = strtrim(getenv('HIRAX_SINGLE_LEAF_ORDERS'));
if ~isempty(order_override)
    tokens = regexp(order_override,'[,\s]+','split');
    settings.surface_orders = str2double(tokens);
end
assert(~isempty(settings.surface_orders) && ...
    all(isfinite(settings.surface_orders)) && ...
    all(settings.surface_orders>=2) && ...
    all(settings.surface_orders==round(settings.surface_orders)) && ...
    all(diff(settings.surface_orders)>0), ...
    ['Surface orders must be a strictly increasing list of ' ...
    'integers >= 2.'])
assert(all(diff(settings.integration_orders)>0) && ...
    settings.integration_orders(end)>max(settings.surface_orders), ...
    'The final integration order must exceed every surface order.')

output_dir = fullfile(test_dir,'results','single_leaf_tip');
if ~isfolder(output_dir), mkdir(output_dir); end
log_file = fullfile(output_dir,'run_single_leaf_tip.log');
diary(log_file)
diary on
diary_cleanup = onCleanup(@()diary('off'));

%% Fix the source independently of surface quadrature order
% Build the four-leaf locator only to obtain the exact outlines and their
% placements. The solve below contains just settings.leaf_id.
locator_options = settings.geometry_options;
locator_options.norder = min(settings.surface_orders);
if isfield(locator_options,'leaf_id')
    locator_options = rmfield(locator_options,'leaf_id');
end
[locator_surface,locator_parts] = generate_leaves_surfer( ...
    settings.thickness,locator_options);
[source_leaf_id,source_point,target_nearest_point,nearest_distance] = ...
    neighboring_source(locator_parts.outlines,locator_parts.rotation, ...
    locator_parts.shift,settings.leaf_id,settings.outline_samples_per_panel);

assert(source_leaf_id ~= settings.leaf_id)
assert(nearest_distance>10*eps(settings.leaf_radius), ...
    'Selected and source leaf outlines touch.')
assert(abs(source_point(3))<100*eps(settings.leaf_radius), ...
    'The source locator must lie on the side-wall equator.')
settings.source_leaf_id = source_leaf_id;
settings.source_point = source_point;
settings.target_nearest_point = target_nearest_point;
settings.nearest_leaf_separation = nearest_distance;
settings.source_info.r = source_point;
settings.source_info.edips = settings.electric_dipole;
settings.locator_npatches = locator_surface.npatches;
clear locator_surface locator_parts locator_options

% These fixed, smooth-field probes are shared with the other leaf
% convergence diagnostics. They stay far enough from the surface for
% independent direct quadrature to be reliable.
settings.probe_targets = settings.leaf_radius*[ ...
    0,0,2,0,-1.6,1.2; ...
    0,0,0,2,.7,-1.3; ...
    2,-2,.4,-.3,.8,-1.1];
settings.probe_labels = {'+z','-z','+x','+y','probe A','probe B'};

fprintf('Single distinct-tip leaf self-convergence\n')
fprintf('geometry revision: %s\n',settings.geometry_revision)
fprintf('target leaf %d; source lies on nearest point of leaf %d\n', ...
    settings.leaf_id,settings.source_leaf_id)
fprintf('source point:         [%.15g %.15g %.15g] mm\n',source_point)
fprintf('nearest target point: [%.15g %.15g %.15g] mm\n',target_nearest_point)
fprintf('outline separation: %.15g mm\n',nearest_distance)
fprintf('surface orders: %s\n',mat2str(settings.surface_orders))

if settings.geometry_only
    save(fullfile(output_dir,'geometry_only.mat'),'settings')
    fprintf('Geometry-only check complete; solver was not run.\n')
    return
end

%% Solve each order and evaluate at identical physical targets
number_of_orders = numel(settings.surface_orders);
levels = cell(1,number_of_orders);
expected_npatches = [];

for order_id = 1:number_of_orders
    norder = settings.surface_orders(order_id);
    level_file = fullfile(output_dir,sprintf('order%02d.mat',norder));
    fprintf('\n========== surface order %d ==========\n',norder)
    try
        geometry_timer = tic;
        opts = settings.geometry_options;
        opts.norder = norder;
        opts.leaf_id = settings.leaf_id;
        [S,parts] = generate_leaves_surfer(settings.thickness,opts);
        geometry_time = toc(geometry_timer);
        assert(strcmp(parts.geometry_revision,settings.geometry_revision), ...
            'Geometry revision changed during the convergence run.')
        if isempty(expected_npatches)
            expected_npatches = S.npatches;
        else
            assert(S.npatches==expected_npatches, ...
                'Patch topology changed between surface orders.')
        end
        sampled_source_distance = min(vecnorm(S.r-source_point,2,1));
        assert(sampled_source_distance >= nearest_distance*(1-1e-7), ...
            'A target-leaf node is closer than the outline separation.')

        % Store only one row from each z-reflected target pair. This uses
        % operator symmetry and does not impose symmetry on the solution.
        reflection = whole_leaf_midplane_permutation( ...
            S,parts,settings.leaf_id);
        correction_target_ids = find((1:S.npts)<reflection);

        fprintf(['Order %d: assembling batched self corrections using ' ...
            'midplane reuse\n'],norder)
        quadrature_timer = tic;
        corrections = nrccie_quad_corr_block(S,settings.eps_quad, ...
            settings.zk,[],settings.quad_batch_size,correction_target_ids);
        quadrature_time = toc(quadrature_timer);
        assert(2*numel(correction_target_ids)==S.npts && ...
            all(cellfun(@(B)size(B,1)==numel(correction_target_ids), ...
            corrections)), ...
            'Reflected correction rows do not match the node pairing.')
        correction_nonzeros = sum(cellfun(@nnz,corrections));
        correction_bytes = cell_storage_bytes(corrections);
        fprintf(['Order %d: %d stored / %d reflected correction nonzeros, ' ...
            '%.2f GiB\n'],norder,correction_nonzeros, ...
            2*correction_nonzeros,correction_bytes/2^30)

        operator.npts = S.npts;
        operator.r = S.r;
        operator.wts = S.wts(:).';
        operator.n = S.n;
        operator.ru = S.du./vecnorm(S.du,2,1);
        operator.rv = cross(operator.n,operator.ru,1);
        operator.zk = settings.zk;
        operator.alpha = settings.alpha;
        operator.eps_fmm = settings.eps_fmm;
        operator.apply_corrections = ...
            @(density)apply_single_corrections( ...
            corrections,density,reflection);
        matvec = @(x)apply_nrccie(x,operator);

        [einc,hinc] = em3d.incoming_sources( ...
            settings.zk,settings.source_info,S,'electric dipole');
        assert(all(isfinite([einc(:);hinc(:)])), ...
            'Order %d produced a nonfinite incident field.',norder)
        normal_einc = sum(S.n.*einc,1);
        tangent_rhs = cross(S.n,hinc,1)-settings.alpha*( ...
            S.n.*normal_einc-einc);
        rhs_components = [sum(operator.ru.*tangent_rhs,1); ...
            sum(operator.rv.*tangent_rhs,1);normal_einc];
        rhs = rhs_components(:);

        fprintf('Order %d: starting GMRES for %d unknowns\n',norder,numel(rhs))
        solve_timer = tic;
        [solution,flag,relres,iter,resvec] = gmres_with_progress( ...
            matvec,rhs,settings.gmres_restart,settings.eps_gmres, ...
            settings.gmres_maxit);
        solve_time = toc(solve_timer);
        true_relative_residual = norm(matvec(solution)-rhs)/norm(rhs);
        assert(flag==0 && true_relative_residual<=5*settings.eps_gmres, ...
            'Order %d GMRES residual %.3e exceeds the requested tolerance.', ...
            norder,true_relative_residual)

        components = reshape(solution,3,S.npts);
        surface_current = operator.ru.*components(1,:)+ ...
            operator.rv.*components(2,:);
        surface_charge = components(3,:);
        density = [surface_current;surface_charge];
        clear matvec operator corrections solution components rhs einc hinc ...
            reflection correction_target_ids normal_einc tangent_rhs

        electric_fields = cell(1,numel(settings.integration_orders));
        integration_times = zeros(1,numel(settings.integration_orders));
        for integration_id = 1:numel(settings.integration_orders)
            integration_order = settings.integration_orders(integration_id);
            integration_timer = tic;
            T = oversample(S,integration_order);
            resampled_density = resample_values(S,density,T);
            electric_fields{integration_id} = direct_electric( ...
                T,resampled_density,settings.probe_targets,settings.zk);
            integration_times(integration_id) = toc(integration_timer);
            assert(all(isfinite(electric_fields{integration_id}(:))), ...
                'Order %d produced a nonfinite probe field.',norder)
            fprintf(['Order %d: integration order %d, %d nodes, ' ...
                'elapsed %.1f s\n'],norder,integration_order,T.npts, ...
                integration_times(integration_id))
            clear T resampled_density
        end
        integration_relative_error = relative_max_difference( ...
            electric_fields{1},electric_fields{end});

        level = struct( ...
            'finished',true,'surface_order',norder, ...
            'npatches',S.npatches,'npts',S.npts, ...
            'sampled_source_distance',sampled_source_distance, ...
            'gmres_flag',flag,'gmres_relative_residual',relres, ...
            'gmres_iterations',iter, ...
            'number_of_iterations',numel(resvec)-1, ...
            'gmres_residual_history',resvec, ...
            'true_relative_residual',true_relative_residual, ...
            'geometry_time',geometry_time, ...
            'quadrature_time',quadrature_time, ...
            'solve_time',solve_time, ...
            'correction_nonzeros',correction_nonzeros, ...
            'reflected_correction_nonzeros',2*correction_nonzeros, ...
            'correction_bytes',correction_bytes, ...
            'midplane_reused',true, ...
            'integration_times',integration_times, ...
            'integration_relative_error',integration_relative_error, ...
            'scattered_electric',electric_fields{end}, ...
            'scattered_electric_by_integration_order',{electric_fields}, ...
            'resolved_geometry_options',parts.options);
        levels{order_id} = level;
        save(level_file,'settings','level','S','surface_current', ...
            'surface_charge','rhs_components','-v7.3')
        fprintf(['order %d: %d nodes, %d iterations, true residual %.3e, ' ...
            'integration change %.3e, solve %.1f s\n'],norder,S.npts, ...
            level.number_of_iterations,true_relative_residual, ...
            integration_relative_error,solve_time)
        clear S parts surface_current surface_charge density
    catch run_error
        failed_level = struct('finished',false,'surface_order',norder, ...
            'error_identifier',run_error.identifier, ...
            'error_message',run_error.message, ...
            'error_report',getReport(run_error,'extended','hyperlinks','off'));
        save(level_file,'settings','failed_level','-v7.3')
        fprintf(2,'order %d failed; saved %s\n%s\n', ...
            norder,level_file,failed_level.error_report)
        rethrow(run_error)
    end
end

%% Self-convergence against the next level and the highest-order reference
orders = settings.surface_orders(:);
npts = cellfun(@(x)x.npts,levels).';
iterations = cellfun(@(x)x.number_of_iterations,levels).';
true_residual = cellfun(@(x)x.true_relative_residual,levels).';
integration_error = cellfun(@(x)x.integration_relative_error,levels).';
successive_error = nan(number_of_orders,1);
reference_error = nan(number_of_orders,1);
reference_digits = nan(number_of_orders,1);

if number_of_orders>1
    reference_field = levels{end}.scattered_electric;
    reference_scale = max(abs(reference_field),[],'all');
    assert(reference_scale>0,'The highest-order probe field is identically zero.')
    for j = 1:number_of_orders-1
        successive_error(j) = relative_max_difference( ...
            levels{j}.scattered_electric,levels{j+1}.scattered_electric);
        reference_error(j) = max(abs( ...
            levels{j}.scattered_electric-reference_field),[],'all')/ ...
            reference_scale;
        reference_digits(j) = -log10(max(reference_error(j),realmin));
    end
else
    reference_scale = max(abs(levels{1}.scattered_electric),[],'all');
    assert(reference_scale>0,'The probe field is identically zero.')
end

convergence_table = table(orders,npts,iterations,true_residual, ...
    integration_error,successive_error,reference_error,reference_digits);
disp(convergence_table)

if number_of_orders>=3
    from_order = orders(1:end-2);
    to_order = orders(2:end-1);
    numerator = max(reference_error(1:end-2),realmin);
    denominator = max(reference_error(2:end-1),realmin);
    error_reduction = numerator./denominator;
    digits_gained_per_order = log10(error_reduction)./(to_order-from_order);
    factor_per_order = 10.^digits_gained_per_order;
else
    from_order = zeros(0,1);
    to_order = zeros(0,1);
    error_reduction = zeros(0,1);
    digits_gained_per_order = zeros(0,1);
    factor_per_order = zeros(0,1);
end
rate_table = table(from_order,to_order,error_reduction, ...
    digits_gained_per_order,factor_per_order);
disp(rate_table)

result.settings = settings;
result.levels = levels;
result.convergence_table = convergence_table;
result.rate_table = rate_table;
result.reference_order = orders(end);
result.reference_scale = reference_scale;
save(fullfile(output_dir,'self_convergence.mat'),'result','-v7.3')
writetable(convergence_table,fullfile(output_dir,'self_convergence.csv'))
writetable(rate_table,fullfile(output_dir,'self_convergence_rates.csv'))
fprintf('Saved convergence results in %s\n',output_dir)


function [source_leaf,source_point,target_point,distance] = ...
    neighboring_source(outlines,rotation,shift,target_leaf,points_per_panel)
[target_samples,target_panel,target_parameter] = sample_outline( ...
    outlines{target_leaf},rotation(:,:,target_leaf),shift(:,target_leaf), ...
    points_per_panel);
source_leaf = NaN;
source_point = nan(3,1);
target_point = nan(3,1);
distance = inf;
for leaf = setdiff(1:4,target_leaf)
    [samples,panel,parameter] = sample_outline( ...
        outlines{leaf},rotation(:,:,leaf),shift(:,leaf),points_per_panel);
    [sample_distance,target_id,source_id] = nearest_pair( ...
        target_samples(1:2,:),samples(1:2,:));
    [candidate_target,candidate_source,candidate_distance] = refine_pair( ...
        outlines{target_leaf},outlines{leaf},rotation(:,:,target_leaf), ...
        rotation(:,:,leaf),shift(:,target_leaf),shift(:,leaf), ...
        target_panel(target_id),panel(source_id), ...
        target_parameter(target_id),parameter(source_id));
    if isfinite(candidate_distance) && candidate_distance<distance
        distance = candidate_distance;
        source_leaf = leaf;
        source_point = candidate_source;
        target_point = candidate_target;
    elseif sample_distance<distance
        distance = sample_distance;
        source_leaf = leaf;
        source_point = samples(:,source_id);
        target_point = target_samples(:,target_id);
    end
end
assert(isfinite(distance),'Could not locate a neighboring source leaf.')
end


function [points,panel_ids,parameters] = sample_outline(outline,R,shift,n)
t = linspace(-1,1,n);
P = reshape(lege.pols(t,outline.order-1),outline.order,[]);
np = outline.number_of_panels;
points = zeros(3,n*np);
panel_ids = zeros(1,n*np);
parameters = repmat(t,1,np);
for panel = 1:np
    ids = (panel-1)*n+(1:n);
    local = outline.position_coefficients(:,:,panel)*P;
    points(:,ids) = R*[local;zeros(1,n)]+shift;
    panel_ids(ids) = panel;
end
end


function [distance,ia,ib] = nearest_pair(a,b)
block_size = 512;
distance_squared = inf;
ia = NaN;
ib = NaN;
aa = sum(a.^2,1).';
for first = 1:block_size:size(b,2)
    ids = first:min(first+block_size-1,size(b,2));
    bb = b(:,ids);
    d2 = max(0,aa+sum(bb.^2,1)-2*(a.'*bb));
    [value,index] = min(d2(:));
    if value<distance_squared
        [ia,j] = ind2sub(size(d2),index);
        ib = ids(j);
        distance_squared = value;
    end
end
distance = sqrt(distance_squared);
end


function [a,b,distance] = refine_pair(oa,ob,Ra,Rb,sa,sb,pa,pb,ta,tb)
objective = @(u)pair_objective(u,oa,ob,Ra,Rb,sa,sb,pa,pb);
options = optimset('Display','off','TolX',1e-13,'TolFun',1e-22, ...
    'MaxFunEvals',800,'MaxIter',400);
u = fminsearch(objective,[ta,tb],options);
u = max(-1,min(1,u));
a = panel_point(oa,pa,u(1),Ra,sa);
b = panel_point(ob,pb,u(2),Rb,sb);
distance = norm(a-b);
end


function value = pair_objective(u,oa,ob,Ra,Rb,sa,sb,pa,pb)
bounded = max(-1,min(1,u));
a = panel_point(oa,pa,bounded(1),Ra,sa);
b = panel_point(ob,pb,bounded(2),Rb,sb);
value = sum((a-b).^2)+1e6*sum((u-bounded).^2);
end


function point = panel_point(outline,panel,t,R,shift)
P = reshape(lege.pols(t,outline.order-1),outline.order,[]);
local = outline.position_coefficients(:,:,panel)*P;
point = R*[local;0]+shift;
end


function permutation = whole_leaf_midplane_permutation(S,parts,leaf_id)
common = parts.common_groups;
local = parts.local_groups{leaf_id};
offset = parts.local_common_surfer.npatches;
groups.top_core = [common.top_core,offset+local.top_core];
groups.top_collar = [common.top_collar,offset+local.top_collar];
groups.upper_wall_by_profile = [common.upper_wall_by_profile, ...
    offset+local.upper_wall_by_profile];
groups.upper_wall = reshape(groups.upper_wall_by_profile,1,[]);
groups.lower_wall_by_profile = [common.lower_wall_by_profile, ...
    offset+local.lower_wall_by_profile];
groups.lower_wall = reshape(groups.lower_wall_by_profile,1,[]);
groups.bottom_collar = [common.bottom_collar,offset+local.bottom_collar];
groups.bottom_core = [common.bottom_core,offset+local.bottom_core];
permutation = leaf_midplane_permutation(S,groups);
end


function [potential,gx,gy,gz] = ...
    apply_single_corrections(B,density,reflection)
values = cell(1,4);
[values{:}] = apply_midplane_quad_corr( ...
    B,density.',reflection,reflection);
potential = values{1}.';
gx = values{2}.';
gy = values{3}.';
gz = values{4}.';
end


function bytes = cell_storage_bytes(values)
bytes = 0;
for k = 1:numel(values)
    value = values{k}; %#ok<NASGU>
    metadata = whos('value');
    bytes = bytes+metadata.bytes;
end
end


function out = resample_values(S,f,T)
% Apply patch interpolation without constructing a global interpolation matrix.
coefficients = vals2coefs(S,f);
out = complex(zeros(size(f,1),T.npts));
[kinds,~,groups] = unique([S.norders(:),S.iptype(:)],'rows');
for group = 1:size(kinds,1)
    patches = find(groups==group);
    norder = kinds(group,1);
    uv = T.uvs_targ(:,T.ixyzs(patches(1)):T.ixyzs(patches(1)+1)-1);
    if kinds(group,2)==1
        B = koorn.pols(norder,uv);
    elseif kinds(group,2)==11
        B = polytens.lege.pols(norder,uv);
    else
        B = polytens.cheb.pols(norder,uv);
    end
    for patch = patches(:).'
        old = S.ixyzs(patch):S.ixyzs(patch+1)-1;
        new = T.ixyzs(patch):T.ixyzs(patch+1)-1;
        out(:,new) = coefficients(:,old)*B;
    end
end
end


function E = direct_electric(S,density,targets,zk)
% Smooth-target NRCCIE electric field: i*k*S[J] - grad(S[rho]).
E = complex(zeros(3,size(targets,2)));
w = S.wts(:).';
for target = 1:size(targets,2)
    delta = targets(:,target)-S.r;
    distance = vecnorm(delta);
    assert(all(distance>0),'A field probe lies on a quadrature node.')
    green = exp(1i*zk*distance)./(4*pi*distance);
    gradient = delta.*(green.*(1i*zk*distance-1)./distance.^2);
    E(:,target) = 1i*zk*sum(density(1:3,:).*(green.*w),2)- ...
        sum(gradient.*(density(4,:).*w),2);
end
end


function error = relative_max_difference(a,b)
scale = max(abs(b),[],'all');
assert(scale>0,'Cannot form a relative error against a zero field.')
error = max(abs(a-b),[],'all')/scale;
end
