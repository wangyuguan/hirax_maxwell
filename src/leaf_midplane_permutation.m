function [permutation,upper_nodes] = leaf_midplane_permutation(S,g)
% Node partners under (x,y,z) -> (x,y,-z), using the leaf patch groups.
upper=[g.top_core,g.top_collar,g.upper_wall];
lower=[g.bottom_core,g.bottom_collar, ...
    reshape(flipud(g.lower_wall_by_profile),1,[])];
assert(numel(upper)==numel(lower) && 2*numel(upper)==S.npatches)
mode=[ones(1,numel(g.top_core)+numel(g.top_collar)), ...
    2*ones(1,numel(g.upper_wall))];
permutation=zeros(1,S.npts);
upper_nodes=zeros(1,S.npts/2);
cache=struct(); count=0;
for j=1:numel(upper)
    a=upper(j); b=lower(j);
    assert(S.norders(a)==S.norders(b) && S.iptype(a)==S.iptype(b))
    ia=S.ixyzs(a):S.ixyzs(a+1)-1;
    ib=S.ixyzs(b):S.ixyzs(b+1)-1;
    key=sprintf('k%d_%d_%d',S.norders(a),S.iptype(a),mode(j));
    if ~isfield(cache,key)
        uv=S.uvs_targ(:,ia);
        if mode(j)==1
            reflected=uv([2,1],:);
        else
            reflected=[-uv(1,:);uv(2,:)];
        end
        distance=(uv(1,:).'-reflected(1,:)).^2+ ...
            (uv(2,:).'-reflected(2,:)).^2;
        [error,p]=min(distance,[],1);
        assert(max(error)<1e-24 && isequal(p(p),1:numel(p)), ...
            'Reference nodes do not have the required reflection symmetry.')
        cache.(key)=p;
    end
    ib=ib(cache.(key));
    permutation(ia)=ib;
    permutation(ib)=ia;
    upper_nodes(count+(1:numel(ia)))=ia;
    count=count+numel(ia);
end
assert(count==S.npts/2 && all(permutation>0) && ...
    isequal(permutation(permutation),1:S.npts))
reflection=[1;1;-1];
scale=max(1,max(abs(S.r),[],'all'));
assert(max(abs(S.r(:,permutation)-reflection.*S.r),[],'all')<2e4*eps(scale), ...
    'Leaf positions are not symmetric about z=0.')
assert(max(abs(S.n(:,permutation)-reflection.*S.n),[],'all')<2e4*eps, ...
    'Leaf normals do not have midplane reflection symmetry.')
weights=S.wts(:).';
assert(max(abs(weights(permutation)-weights))<2e4*eps(max(abs(weights))), ...
    'Leaf quadrature weights do not have midplane reflection symmetry.')
end
