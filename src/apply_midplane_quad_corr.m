function [p,x,y,z] = apply_midplane_quad_corr(B,d,source_perm,target_perm)
% B contains one target from each z-reflected pair; d is nsource-by-ndensity.
rows=find((1:numel(target_perm))<target_perm);
nd=size(d,2);
d=[d,d(source_perm,:)];
out=cell(1,4); signs=[1,1,1,-1];
for k=1:4
    v=B{k}*d;
    out{k}=complex(zeros(numel(target_perm),nd));
    out{k}(rows,:)=v(:,1:nd);
    out{k}(target_perm(rows),:)=signs(k)*v(:,nd+(1:nd));
end
[p,x,y,z]=out{:};
end
