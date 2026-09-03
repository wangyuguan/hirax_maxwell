function [Cslp,Cx,Cy,Cz,info] = fourfold_quad_corr_mats( ...
    S1,S2,S3,S4,eps_quad,zk)
%FOURFOLD_QUAD_CORR_MATS Reuse one block row for four rotated surfaces.
%
% The surfaces must form the clockwise rotation cycle S1 -> S2 -> S4 ->
% S3 -> S1. The returned matrices use the merge order [S1,S2,S3,S4].
% Only geometry/operator symmetry is used; densities on the four surfaces
% may be arbitrary and unrelated.

surfaces = {S1,S2,S4,S3};
cycle_to_merge = [1,2,4,3];
Q = [0,1,0;-1,0,0;0,0,1];

npts = S1.npts;
for isurface = 1:4
    assert(surfaces{isurface}.npts == npts, ...
        'All four surfaces must have the same number of nodes.')
end

scale = max(1,max(vecnorm(S1.r,2,1)));
symmetry_tolerance = 1e3*eps(scale);
position_errors = zeros(1,4);
normal_errors = zeros(1,4);
weight_errors = zeros(1,4);
for isurface = 1:4
    inext = mod(isurface,4)+1;
    position_errors(isurface) = max(vecnorm( ...
        Q*surfaces{isurface}.r-surfaces{inext}.r,2,1));
    normal_errors(isurface) = max(vecnorm( ...
        Q*surfaces{isurface}.n-surfaces{inext}.n,2,1));
    weight_errors(isurface) = max(abs( ...
        surfaces{isurface}.wts-surfaces{inext}.wts));
end
assert(max(position_errors) < symmetry_tolerance && ...
    max(normal_errors) < symmetry_tolerance && ...
    max(weight_errors) < symmetry_tolerance, ...
    'The four surfaces do not have the required rotational symmetry.')

% Construct one target block row. The self block uses singular
% self-to-self quadrature; the other blocks are close-surface corrections.
Bslp = cell(1,4);
Bx = cell(1,4);
By = cell(1,4);
Bz = cell(1,4);
for source_cycle = 1:4
    source_surface = surfaces{source_cycle};
    if source_cycle == 1
        Bslp{source_cycle} = em3d.slp.get_quad_corr_mat( ...
            source_surface,eps_quad,zk);
        [Bx{source_cycle},By{source_cycle},Bz{source_cycle}] = ...
            em3d.sgrad.get_quad_corr_mat(source_surface,eps_quad,zk);
    else
        Bslp{source_cycle} = em3d.slp.get_quad_corr_mat( ...
            source_surface,eps_quad,zk,S1);
        [Bx{source_cycle},By{source_cycle},Bz{source_cycle}] = ...
            em3d.sgrad.get_quad_corr_mat( ...
            source_surface,eps_quad,zk,S1);
    end
end

% Rotate and copy the block row into the full merge ordering. The scalar
% correction is invariant; its target derivatives rotate as a vector.
Cslp_blocks = cell(4,4);
Cx_blocks = cell(4,4);
Cy_blocks = cell(4,4);
Cz_blocks = cell(4,4);
for target_cycle = 1:4
    Qtarget = Q^(target_cycle-1);
    target_merge = cycle_to_merge(target_cycle);
    for source_cycle = 1:4
        relative_cycle = mod(source_cycle-target_cycle,4)+1;
        source_merge = cycle_to_merge(source_cycle);

        Cslp_blocks{target_merge,source_merge} = Bslp{relative_cycle};
        Cx_blocks{target_merge,source_merge} = ...
            Qtarget(1,1)*Bx{relative_cycle}+ ...
            Qtarget(1,2)*By{relative_cycle}+ ...
            Qtarget(1,3)*Bz{relative_cycle};
        Cy_blocks{target_merge,source_merge} = ...
            Qtarget(2,1)*Bx{relative_cycle}+ ...
            Qtarget(2,2)*By{relative_cycle}+ ...
            Qtarget(2,3)*Bz{relative_cycle};
        Cz_blocks{target_merge,source_merge} = ...
            Qtarget(3,1)*Bx{relative_cycle}+ ...
            Qtarget(3,2)*By{relative_cycle}+ ...
            Qtarget(3,3)*Bz{relative_cycle};
    end
end

Cslp = cell2mat(Cslp_blocks);
Cx = cell2mat(Cx_blocks);
Cy = cell2mat(Cy_blocks);
Cz = cell2mat(Cz_blocks);

info = struct();
info.cycle_to_merge = cycle_to_merge;
info.rotation = Q;
info.position_errors = position_errors;
info.normal_errors = normal_errors;
info.weight_errors = weight_errors;
info.base_block_nonzeros = zeros(4,4);
for source_cycle = 1:4
    info.base_block_nonzeros(source_cycle,:) = [ ...
        nnz(Bslp{source_cycle}),nnz(Bx{source_cycle}), ...
        nnz(By{source_cycle}),nnz(Bz{source_cycle})];
end
end
