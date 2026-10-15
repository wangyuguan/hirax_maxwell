function diagnose_single_leaf
% Diagnose existing results without solving again or creating figures.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
orders = [4 6 8 10 12];
results = struct([]);
previous = [];
for k = 1:numel(orders)
    clock0 = tic;
    d = load(fullfile(root,'data',sprintf( ...
        'run_hirax_single_leaf_order%d.mat',orders(k))));
    S = d.S;
    c = d.case_settings;
    source = struct('r',c.x0,'edips',-c.p0);
    [E,H] = em3d.incoming_sources(c.zk,source,S,'electric dipole');
    rhs = make_rhs(S,E,H,c.alpha);
    a = struct('order',orders(k));
    a.incident_saved_difference = max(abs([E-d.einc;H-d.hinc]),[],'all');
    a.area = sum(S.wts);
    a.volume = sum(sum(S.r.*S.n,1).*S.wts(:).')/3;
    a.gmres_residual = d.level.relative_residual;
    a.gmres_iterations = d.level.iterations;
    a.position_mismatch = d.level.max_position_mismatch;
    a.normal_mismatch_degrees = d.level.max_normal_mismatch_degrees;
    a.J_L2 = weighted_norm(d.densities(1:3,:),S.wts);
    a.rho_L2 = weighted_norm(d.densities(4,:),S.wts);
    a.normal_current_relative = max(abs(sum(S.n.*d.densities(1:3,:),1))) / ...
        max(abs(d.densities(1:3,:)),[],'all');
    rhs_scale = [max(abs(rhs(1:2,:)),[],'all'); ...
        max(abs(rhs(1:2,:)),[],'all');max(abs(rhs(3,:)))];
    a.rhs_tail = max(surf_fun_error(S,rhs./rhs_scale),[],2).';
    density_scale = max(abs(d.densities),[],2);
    a.density_tail = max(surf_fun_error(S,d.densities./density_scale),[],2).';
    a.saved_E = d.level.escat;
    a.direct_native_E = direct_electric(S,d.densities,c.test_targets,c.zk);
    fprintf('\nORDER %d: geometry/coefficients loaded; checking independent nodes\n',orders(k));
    So = oversample(S,16);
    fo = resample_values(S,[rhs;d.densities],So);
    [Eo,Ho] = em3d.incoming_sources(c.zk,source,So,'electric dipole');
    exact_rhs = make_rhs(So,Eo,Ho,c.alpha);
    rhs_error = fo(1:3,:)-exact_rhs;
    a.rhs_check_max_relative = max(abs(rhs_error)./rhs_scale,[],2).';
    a.rhs_check_L2_relative = weighted_norm(rhs_error,So.wts)/ ...
        weighted_norm(exact_rhs,So.wts);
    [~,worst_node] = max(max(abs(rhs_error)./rhs_scale,[],1));
    a.rhs_worst_patch = So.patch_id(worst_node);
    a.rhs_worst_position = So.r(:,worst_node).';
    a.direct_order16_E = direct_electric(So,fo(4:7,:),c.test_targets,c.zk);
    a.saved_vs_direct_relative = max(abs(a.saved_E-a.direct_order16_E),[],'all')/ ...
        max(abs(a.direct_order16_E),[],'all');
    a.native_vs_order16_relative = max(abs(a.direct_native_E-a.direct_order16_E),[],'all')/ ...
        max(abs(a.direct_order16_E),[],'all');
    a.previous_geometry_difference = NaN;
    a.previous_density_L2_relative = NaN;
    a.previous_direct_E_relative = NaN;
    if ~isempty(previous)
        a.previous_geometry_difference = max(vecnorm(So.r-previous.r));
        a.previous_density_L2_relative = weighted_norm(fo(4:7,:)-previous.densities,So.wts)/ ...
            weighted_norm(fo(4:7,:),So.wts);
        a.previous_direct_E_relative = max(abs(a.direct_order16_E-previous.E),[],'all')/ ...
            max(abs(a.direct_order16_E),[],'all');
    end
    previous = struct('r',So.r,'densities',fo(4:7,:),'E',a.direct_order16_E);
    a.elapsed = toc(clock0);
    if k==1
        results = a;
    else
        results(k) = a;
    end
    disp(a)
    save(fullfile(root,'diagnostics','single_leaf_diagnosis.mat'),'results');
end
fprintf('DIAGNOSIS_DONE\n');
end

function rhs = make_rhs(S,E,H,alpha)
u = S.du./vecnorm(S.du);
v = cross(S.n,u,1);
nE = sum(S.n.*E,1);
t = cross(S.n,H,1)-alpha*(S.n.*nE-E);
rhs = [sum(u.*t,1);sum(v.*t,1);nE];
end

function val = weighted_norm(f,w)
val = sqrt(sum(sum(abs(f).^2,1).*w(:).'));
end

function out = resample_values(S,f,T)
% Apply patch interpolation in batches without constructing a large matrix.
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
% Independent smooth quadrature, valid here because probes are far away.
% E = ik S_k[J] - grad S_k[rho], G_k = exp(ik r)/(4 pi r).
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
