% Compare the four scalar kernels with the native nine-kernel NRCCIE.
root = fileparts(fileparts(mfilename('fullpath')));
run('../fmm3dbie-hirax-dev/matlab/startup.m')
addpath(genpath('src'))
addpath('../fmm3dbie-hirax-dev/FMM3D/matlab')
clear pth dir
orders = [4,8,12]; zk = .7; alpha = .8; eps_quad = 1e-10;
charge_errors = zeros(size(orders)); current_errors = charge_errors;
cube_errors = zeros(1,2); rng(8)
for io = 0:numel(orders)
    if io==0
        S = flat_cube(4);
    else
        S = geometries.sphere(1,1,[0;0;0],orders(io),11);
    end
    N = S.npts; P = S.npatches;
    [srcvals,srccoefs,norders,ixyzs,iptype] = extract_arrays(S);
    % Include every source patch so the reference has no FMM or far rule.
    row_ptr = (1:P:P*N+1).';
    col_ind = repmat((1:P).',N,1);
    iquad = [1;1+cumsum(repmat(diff(ixyzs(:)),N,1))];
    nquad = iquad(end)-1; nnear = numel(col_ind);
    wnear = complex(zeros(9,nquad));
    mex_id_ = 'getnearquad_em_nrccie_pec(i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[xx], i double[xx], i double[x], i dcomplex[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i int64_t[x], i double[x], i int64_t[x], io dcomplex[xx])';
    wnear = fmm3dbie_routs(mex_id_,P,norders,ixyzs,iptype,N,srccoefs,srcvals, ...
        eps_quad,complex([zk;alpha]),1,nnear,row_ptr,col_ind,iquad,1.25,nquad,wnear, ...
        1,P,P+1,P,1,9,N,12,N,1,2,1,1,N+1,nnear,nnear+1,1,1,9,nquad);
    K = cell(3,3);
    for a = 1:3
        for b = 1:3
            K{a,b} = reshape(wnear(3*(a-1)+b,:),N,N).';
        end
    end
    clear wnear
    % Force the scalar correction to cover the same complete patch list.
    opts = struct('rfac',10,'rfac0',1.25);
    C = cell(1,4);
    C{1} = em3d.slp.get_quad_corr_mat(S,eps_quad,zk,opts);
    [C{2:4}] = em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk,opts);
    op = struct('npts',N,'r',S.r,'wts',S.wts(:).','n',S.n, ...
        'ru',S.du./vecnorm(S.du,2,1),'zk',zk,'alpha',alpha,'eps_fmm',1e-12);
    op.rv = cross(S.n,op.ru,1);
    op.apply_corrections = @(d)full_apply(C,d);
    for test = 1:2
        d = complex(zeros(3,N));
        if test==1
            d(3,:) = 1+.2*S.r(1,:)+.3i*S.r(2,:);
        elseif io==0
            face_values = randn(2,6)+1i*randn(2,6);
            d(1:2,:) = face_values(:,S.patch_id);
        else
            J = cross(S.n,repmat([1;.3i;.2],1,N),1);
            d(1,:) = sum(op.ru.*J,1); d(2,:) = sum(op.rv.*J,1);
        end
        ref = .5*d;
        for a = 1:3
            for b = 1:3
                ref(a,:) = ref(a,:)+(K{a,b}*d(b,:).').';
            end
        end
        got = reshape(apply_nrccie(d(:),op),3,N);
        err = norm(got-ref,'fro')/norm(ref,'fro');
        if io==0
            cube_errors(test) = err;
        elseif test==1
            charge_errors(io) = err;
        else
            current_errors(io) = err;
        end
    end
    if io==0
        fprintf('Flat cube, N %d: native charge %.3e; constant-face current %.3e\n',N,cube_errors)
    else
        fprintf('Sphere order %d, N %d: native charge %.3e; smooth current %.3e\n', ...
            orders(io),N,charge_errors(io),current_errors(io))
    end
end
save(fullfile(root,'diagnostics','nrccie_native_check.mat'), ...
    'orders','charge_errors','current_errors','cube_errors','zk','alpha','eps_quad')
% Cartesian-current and local-basis interpolation differ at finite order.
assert(max(cube_errors)<1e-7,'Native kernels disagree with constant tangent bases.')
assert(max(charge_errors)<1e-5,'Native charge blocks disagree.')
assert(current_errors(end)<current_errors(1)/10 && current_errors(end)<1e-4, ...
    'Current representation discrepancy does not decrease with order.')
disp('NATIVE_NRCCIE_PASS')

function [p,x,y,z] = full_apply(C,d)
p = (C{1}*d.').'; x = (C{2}*d.').';
y = (C{3}*d.').'; z = (C{4}*d.').';
end

function S = flat_cube(order)
uv = polytens.lege.nodes(order); np = size(uv,2);
normals = [eye(3),-eye(3)]; vals = zeros(12,6*np);
for f = 1:6
    n = normals(:,f); u = [1;0;0];
    if abs(n(1))==1, u = [0;1;0]; end
    v = cross(n,u); ii = (f-1)*np+(1:np);
    vals(:,ii) = [n+u*uv(1,:)+v*uv(2,:); ...
        repmat(u,1,np);repmat(v,1,np);repmat(n,1,np)];
end
S = surfer(6,order,vals,11);
end
