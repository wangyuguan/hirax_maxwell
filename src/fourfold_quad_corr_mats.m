function [corrections,info] = fourfold_quad_corr_mats( ...
    S1,S2,S3,S4,eps_quad,zk,batch_size,reflection)
%FOURFOLD_QUAD_CORR_MATS Base quadrature blocks for four rotated surfaces.
%
% The surfaces must form the clockwise rotation cycle S1 -> S2 -> S4 ->
% S3 -> S1. Only geometry/operator symmetry is used; densities on the
% four surfaces may be arbitrary and unrelated.
%
% The full correction matrices are block circulant in the cycle ordering:
% the block coupling target cycle index it to source cycle index is holds
% only on mod(is-it,4). This routine returns that single block row rather
% than the assembled matrices, which is a factor of four less storage.
% Use FOURFOLD_APPLY_QUAD_CORR to apply it, or FOURFOLD_EXPAND_QUAD_CORR
% to recover the assembled matrices in the merge order [S1,S2,S3,S4].
%
% The returned struct has fields slp, gx, gy, gz, each a 1-by-4 cell of
% sparse blocks whose target surface is S1 and whose source surface sits
% at the matching relative position in the cycle.
% Optional REFLECTION maps nodes under z -> -z and halves the stored rows.

if nargin<7, batch_size=[]; end
if nargin<8, reflection=[]; end
target_ids=[];
if ~isempty(reflection)
    target_ids=find((1:S1.npts)<reflection);
    if isempty(batch_size), batch_size=2000; end
end
surfaces = {S1,S2,S4,S3};
cycle_to_merge = [1,2,4,3];
Q = [0,1,0;-1,0,0;0,0,1];

npts = S1.npts;
for isurface = 1:4
    assert(surfaces{isurface}.npts == npts, ...
        'All four surfaces must have the same number of nodes.')
    assert(isequal(surfaces{isurface}.norders,S1.norders) && ...
        isequal(surfaces{isurface}.iptype,S1.iptype), ...
        'Rotational reuse requires matching patch orders and types.')
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
    if ~isempty(batch_size)
        target=S1;
        if source_cycle==1, target=[]; end
        B=nrccie_quad_corr_block(source_surface,eps_quad,zk,target,batch_size,target_ids);
        [Bslp{source_cycle},Bx{source_cycle},By{source_cycle},Bz{source_cycle}]=B{:};
    elseif source_cycle == 1
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

% Surfaces on opposite corners have no near interaction at all, so the
% matching block row entry is structurally empty and its multiply is
% skipped when the correction is applied.
nonempty = false(1,4);
for source_cycle = 1:4
    nonempty(source_cycle) = nnz(Bslp{source_cycle})+ ...
        nnz(Bx{source_cycle})+nnz(By{source_cycle})+ ...
        nnz(Bz{source_cycle}) > 0;
end

corrections = struct();
corrections.slp = Bslp;
corrections.gx = Bx;
corrections.gy = By;
corrections.gz = Bz;
corrections.nonempty = nonempty;
corrections.cycle_to_merge = cycle_to_merge;
corrections.rotation = Q;
corrections.component_npts = npts;
corrections.npts = 4*npts;
corrections.reflection = reflection;

info = struct();
info.cycle_to_merge = cycle_to_merge;
info.rotation = Q;
info.position_errors = position_errors;
info.normal_errors = normal_errors;
info.weight_errors = weight_errors;
info.nonempty = nonempty;
info.base_block_nonzeros = zeros(4,4);
for source_cycle = 1:4
    info.base_block_nonzeros(source_cycle,:) = [ ...
        nnz(Bslp{source_cycle}),nnz(Bx{source_cycle}), ...
        nnz(By{source_cycle}),nnz(Bz{source_cycle})];
end
info.stored_nonzeros = sum(info.base_block_nonzeros,'all');
info.assembled_nonzeros = 4*(1+~isempty(reflection))*info.stored_nonzeros;
end
