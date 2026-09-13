function [Cslp,Cx,Cy,Cz] = fourfold_expand_quad_corr(corrections)
%FOURFOLD_EXPAND_QUAD_CORR Assemble the full four-surface correction.
%
% Rebuilds the sixteen-block matrices in the merge order [S1,S2,S3,S4]
% from the single block row returned by FOURFOLD_QUAD_CORR_MATS. This
% costs four times the storage (eight with midplane reflection) and is for testing
% and debugging; the solver should call FOURFOLD_APPLY_QUAD_CORR instead.

cycle_to_merge = corrections.cycle_to_merge;
Q = corrections.rotation;
if isfield(corrections,'reflection') && ~isempty(corrections.reflection)
    perm=corrections.reflection; n=numel(perm);
    rows=find((1:n)<perm); names={'slp','gx','gy','gz'}; signs=[1,1,1,-1];
    [~,row_order]=sort([rows,perm(rows)]);
    for k=1:4
        for j=1:4
            B=corrections.(names{k}){j};
            A=[B;signs(k)*B(:,perm)];
            corrections.(names{k}){j}=A(row_order,:);
        end
    end
end

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

        Cslp_blocks{target_merge,source_merge} = ...
            corrections.slp{relative_cycle};
        Cx_blocks{target_merge,source_merge} = ...
            Qtarget(1,1)*corrections.gx{relative_cycle}+ ...
            Qtarget(1,2)*corrections.gy{relative_cycle}+ ...
            Qtarget(1,3)*corrections.gz{relative_cycle};
        Cy_blocks{target_merge,source_merge} = ...
            Qtarget(2,1)*corrections.gx{relative_cycle}+ ...
            Qtarget(2,2)*corrections.gy{relative_cycle}+ ...
            Qtarget(2,3)*corrections.gz{relative_cycle};
        Cz_blocks{target_merge,source_merge} = ...
            Qtarget(3,1)*corrections.gx{relative_cycle}+ ...
            Qtarget(3,2)*corrections.gy{relative_cycle}+ ...
            Qtarget(3,3)*corrections.gz{relative_cycle};
    end
end

Cslp = cell2mat(Cslp_blocks);
Cx = cell2mat(Cx_blocks);
Cy = cell2mat(Cy_blocks);
Cz = cell2mat(Cz_blocks);
end
