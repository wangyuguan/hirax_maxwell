function check_multiple_leaves_convergence(orders)
% No-plot probe-field self-convergence of the saved multiple-leaf solutions.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
if nargin==0
    orders = [4 6 8 10];
end
quadrature_orders = [16 20];
results = cell(size(orders));
for k = 1:numel(orders)
    d = load(fullfile(root,'data',sprintf( ...
        'run_multiple_leaves_order%d.mat',orders(k))));
    assert(d.solver_data.gmres_flag==0,'GMRES did not converge');
    assert(d.settings.norder==orders(k));
    settings = rmfield(d.settings,'norder');
    if k==1
        common_settings = settings;
        % Same six probe coordinates as the earlier single-leaf comparison.
        targets = settings.leaf_radius*[ ...
            0 0 2 0 -1.6 1.2;0 0 0 2 0.7 -1.3;2 -2 0.4 -0.3 0.8 -1.1];
    else
        assert(isequaln(rmfield(settings,'shifts'), ...
            rmfield(common_settings,'shifts')),'Saved settings differ');
        shift_difference = max(abs(settings.shifts-common_settings.shifts),[],'all');
        assert(shift_difference<=64*eps(settings.geometry_radius), ...
            'Leaf placements differ beyond floating-point roundoff');
        fprintf('order=%d placement_difference=%.3e\n',orders(k),shift_difference);
    end
    density = [d.surface_current;d.surface_charge];
    assert(isequal(size(density),[4 d.S.npts]));
    assert(all(isfinite(density(:))));
    r = struct('order',orders(k),'npts',d.S.npts, ...
        'gmres_iterations',d.solver_data.number_of_iterations, ...
        'gmres_residual',d.solver_data.true_relative_residual);
    r.probe_node_distances = zeros(1,size(targets,2));
    for t = 1:size(targets,2)
        r.probe_node_distances(t) = min(vecnorm(targets(:,t)-d.S.r));
    end
    r.electric_fields = cell(size(quadrature_orders));
    for q = 1:numel(quadrature_orders)
        timer = tic;
        T = oversample(d.S,quadrature_orders(q));
        f = resample_values(d.S,density,T);
        r.electric_fields{q} = direct_electric(T,f,targets,settings.zk);
        fprintf('order=%d integration_order=%d nodes=%d elapsed=%.2fs\n', ...
            orders(k),quadrature_orders(q),T.npts,toc(timer));
        clear T f
    end
    r.integration_check_relative = max(abs( ...
        r.electric_fields{1}-r.electric_fields{2}),[],'all')/ ...
        max(abs(r.electric_fields{2}),[],'all');
    results{k} = r;
    fprintf('order=%d GMRES=%d residual=%.8e integration_change=%.8e\n', ...
        r.order,r.gmres_iterations,r.gmres_residual,r.integration_check_relative);
end
reference_order = orders(end);
reference = results{end}.electric_fields{end};
reference_max = max(abs(reference),[],'all');
comparison_orders = orders(1:end-1).';
absolute_difference = zeros(size(comparison_orders));
for k = 1:numel(comparison_orders)
    absolute_difference(k) = max(abs( ...
        results{k}.electric_fields{end}-reference),[],'all');
end
relative_difference = absolute_difference/reference_max;
digits = -log10(relative_difference);
digits_per_order = diff(digits)./diff(comparison_orders);
reduction = absolute_difference(1:end-1)./absolute_difference(2:end);
convergence = table(comparison_orders,absolute_difference,relative_difference,digits);
fprintf('REFERENCE order=%d max_field=%.16e\n',reference_order,reference_max);
disp(convergence)
for k = 1:numel(digits_per_order)
    fprintf(['RATE %d->%d ERROR_REDUCTION %.12g ' ...
        'DIGITS_PER_ORDER %.12g FACTOR_PER_ORDER %.12g\n'], ...
        comparison_orders(k),comparison_orders(k+1), ...
        reduction(k),digits_per_order(k),10.^digits_per_order(k));
end
fprintf('PROBE_COORDINATES\n'); disp(targets)
save(fullfile(root,'diagnostics','multiple_leaves_convergence.mat'), ...
    'results','targets','common_settings','quadrature_orders', ...
    'reference_order','reference_max','convergence','digits_per_order','reduction');
fprintf('MULTIPLE_LEAVES_CONVERGENCE_DONE\n');
end

function out = resample_values(S,f,T)
coefficients = vals2coefs(S,f);
out = complex(zeros(size(f,1),T.npts));
[kinds,~,groups] = unique([S.norders(:),S.iptype(:)],'rows');
for g = 1:size(kinds,1)
    patches = find(groups==g);
    p = kinds(g,1);
    uv = T.uvs_targ(:,T.ixyzs(patches(1)):T.ixyzs(patches(1)+1)-1);
    if kinds(g,2)==1
        B = koorn.pols(p,uv);
    elseif kinds(g,2)==11
        B = polytens.lege.pols(p,uv);
    else
        B = polytens.cheb.pols(p,uv);
    end
    for patch = patches(:).'
        old = S.ixyzs(patch):S.ixyzs(patch+1)-1;
        new = T.ixyzs(patch):T.ixyzs(patch+1)-1;
        out(:,new) = coefficients(:,old)*B;
    end
end
end

function E = direct_electric(S,density,targets,zk)
E = complex(zeros(3,size(targets,2)));
w = S.wts(:).';
for t = 1:size(targets,2)
    delta = targets(:,t)-S.r;
    r = vecnorm(delta);
    G = exp(1i*zk*r)./(4*pi*r);
    gradG = delta.*(G.*(1i*zk*r-1)./r.^2);
    E(:,t) = 1i*zk*sum(density(1:3,:).*(G.*w),2)- ...
        sum(gradG.*(density(4,:).*w),2);
end
end
