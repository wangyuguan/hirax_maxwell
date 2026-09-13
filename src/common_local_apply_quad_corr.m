function [potential,gx,gy,gz] = common_local_apply_quad_corr(C,density)
% Apply corrections in the full leaf merge order, without expanding U-U.
assert(size(density,2)==C.npts)
out=cell(1,4);
base=cell(1,4);
[base{:}]=fourfold_apply_quad_corr(C.common,density(:,C.common_nodes));
for k=1:4
    out{k}=complex(zeros(size(density)));
    out{k}(:,C.common_nodes)=base{k};
end
for target=1:8
    ti=C.nodes{target};
    for source=1:8
        B=C.blocks{target,source};
        if isempty(B), continue; end
        d=density(:,C.nodes{source}).';
        if isfield(C,'reflection') && ~isempty(C.reflection)
            v=cell(1,4);
            [v{:}]=apply_midplane_quad_corr(B,d,C.reflection{source},C.reflection{target});
            for k=1:4, out{k}(:,ti)=out{k}(:,ti)+v{k}.'; end
            continue
        end
        for k=1:4
            out{k}(:,ti)=out{k}(:,ti)+(B{k}*d).';
        end
    end
end
[potential,gx,gy,gz]=out{:};
end
