% Proxy-assisted NRCCIE factorization: independent far blocks and self solves.
% The held-out ID check must compress actual DOFs, even if this small mesh
% has too few DOFs for FLAM to eliminate any on its proxy-enabled levels.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'..','FLAM')))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(57)

S = geometries.ellipsoid([2,0.6,0.4],[2,2,2],[0;0;0],7,11);
assert(S.npts==1536)
ndof = 3*S.npts;
tol = 1e-8;
max_nodes = 64;
opts = struct('zk',0.9,'alpha',1,'jump',-0.5,'fmm',false);
coordinates = repelem(S.r,1,3);

% A box surrounding the right end cap. All validation points are outside
% the proxy ellipsoid and none of their rows is included in the kernel ID.
cap_nodes = find(S.r(1,:)>0.9);
slf = reshape(3*cap_nodes+(-2:0).',1,[]);
lower = min(S.r(:,cap_nodes),[],2);
upper = max(S.r(:,cap_nodes),[],2);
box_size = upper-lower;
center = (upper+lower)/2;
far_nodes = find(sum(((S.r-center)./box_size).^2,1)>1.65^2);
far = reshape(3*far_nodes+(-2:0).',1,[]);
assert(~isempty(far) && isempty(intersect(slf,far)))
heldout = [nrccie_system_block(S,[],[],opts,far,slf); ...
    nrccie_system_block(S,[],[],opts,slf,far)'];
proxy_orders = [16,32];
proxy_errors = zeros(2,2);
for pass = 1:2
    proxy_opts = opts;
    proxy_opts.proxy_order = proxy_orders(pass);
    callback = nrccie_rskelf_proxy(S,[],proxy_opts,tol);
    [sampling,neighbors] = callback(coordinates,slf,[],box_size,center);
    assert(isempty(neighbors))
    [sk,rd,interpolation] = id(sampling,tol);
    assert(~isempty(rd),'Held-out proxy test did not compress any DOFs.')
    approximation = heldout(:,sk)*interpolation;
    proxy_errors(pass,1) = relative_error( ...
        approximation(1:numel(far),:),heldout(1:numel(far),rd));
    proxy_errors(pass,2) = relative_error( ...
        approximation(numel(far)+1:end,:),heldout(numel(far)+1:end,rd));
    fprintf(['Proxy order %d: %d/%d skeleton DOFs, %d independent far DOFs; ' ...
        'outgoing %.3e, incoming %.3e.\n'], ...
        proxy_orders(pass),numel(sk),numel(slf),numel(far),proxy_errors(pass,:))
end
assert(max(proxy_errors(2,:))<50*tol)

% A remote correction is not a smooth kernel interaction. Add one in each
% direction, involving a DOF that the smooth proxy initially eliminated.
% Recompress and check the actual modified operator on all held-out nodes.
remote_node = far_nodes(1);
local_node = ceil(slf(rd(end))/3);
remote_correction = repmat({sparse(S.npts,S.npts)},1,4);
remote_correction{1}(remote_node,local_node) = 0.04+0.02i;
remote_correction{2}(local_node,remote_node) = 0.07-0.03i;
callback = nrccie_rskelf_proxy(S,remote_correction,proxy_opts,tol);
[sampling,neighbors] = callback(coordinates,slf,[],box_size,center);
assert(isempty(neighbors),'Correction rows must not add inactive tree neighbors.')
[sk,rd,interpolation] = id(sampling,tol);
assert(~isempty(rd))
corrected_holdout = [ ...
    nrccie_system_block(S,[],remote_correction,opts,far,slf); ...
    nrccie_system_block(S,[],remote_correction,opts,slf,far)'];
remote_error = relative_error( ...
    corrected_holdout(:,sk)*interpolation,corrected_holdout(:,rd));
fprintf('Remote correction held-out error: %.3e.\n',remote_error)
assert(remote_error<50*tol)
clear sampling heldout corrected_holdout callback interpolation

% Corrected self operators, both jump signs, compared with independent FMM
% application and dense LU. Reuse the same correction and random probes.
fprintf('Building quadrature corrections for %d nodes ...\n',S.npts)
C = nrccie_quad_corr_block(S,1e-10,opts.zk,[],128);
[~,Aminus] = nrccie_apply(S,zeros(ndof,1),C,opts);
q = randn(ndof,2)+1i*randn(ndof,2);
q = q./vecnorm(q);
results = struct([]);
jumps = [-0.5,0.5];
for pass = 1:2
    opts.jump = jumps(pass);
    opts.verb = 1;
    if pass==1, A = Aminus;
    else, A = Aminus+speye(ndof); end
    fprintf('Factoring jump=%+.1f with proxy, tol %.1e ...\n',opts.jump,tol)
    timer = tic;
    F = nrccie_rskelf(S,C,opts,max_nodes,tol);
    results(pass).factor_seconds = toc(timer);
    fmm_opts = opts;
    fmm_opts.fmm = true;
    fmm_opts.eps_fmm = 1e-11;
    Aq = complex(zeros(size(q)));
    solution = rskelf_sv(F,q);
    Asolution = complex(zeros(size(q)));
    for column = 1:size(q,2)
        Aq(:,column) = nrccie_apply(S,q(:,column),C,fmm_opts);
        Asolution(:,column) = nrccie_apply(S,solution(:,column),C,fmm_opts);
    end
    results(pass).jump = opts.jump;
    results(pass).fmm_dense_error = relative_error(Aq,A*q);
    results(pass).multiply_error = relative_error(rskelf_mv(F,q),Aq);
    results(pass).solve_error = relative_error(solution,A\q);
    results(pass).solve_residual = relative_error(Asolution,q);
    results(pass).below_root = eliminated_dofs(F,F.lvp(end-1));
    results(pass).proxy_levels = eliminated_dofs(F,F.lvp(max(1,end-2)));
    assert(results(pass).below_root>0, ...
        'Full factorization test did not eliminate any DOFs below the root.')
    assert(results(pass).fmm_dense_error<5e-8)
    assert(max([results(pass).multiply_error,results(pass).solve_error, ...
        results(pass).solve_residual])<50*tol)
    fprintf(['Jump %+.1f: %.2f s; FMM/dense %.3e, mv %.3e, ' ...
        'sv/dense %.3e, true residual %.3e; ' ...
        'eliminated below root %d, on proxy levels %d.\n'], ...
        opts.jump,results(pass).factor_seconds,results(pass).fmm_dense_error, ...
        results(pass).multiply_error,results(pass).solve_error, ...
        results(pass).solve_residual,results(pass).below_root, ...
        results(pass).proxy_levels)
    clear F

    if pass==1
        plain_opts = opts;
        plain_opts.proxy = false;
        timer = tic;
        Fplain = nrccie_rskelf(S,C,plain_opts,max_nodes,tol);
        plain_seconds = toc(timer);
        plain_error = relative_error(rskelf_mv(Fplain,q),Aq);
        assert(plain_error<50*tol)
        fprintf('No-proxy comparison: %.2f s; mv/FMM %.3e.\n', ...
            plain_seconds,plain_error)
        clear Fplain
    end
end
fprintf('NRCCIE_RSKELF_PROXY_PASS\n')


function count = eliminated_dofs(F,last)
count = 0;
for j = 1:last, count = count+numel(F.factors(j).rd); end
end


function error = relative_error(a,b)
error = norm(a(:)-b(:))/max(norm(b(:)),realmin);
end
