function [permutation,upper_nodes] = leaf_midplane_permutation(S,g)
% Node partners under (x,y,z) -> (x,y,-z), using the leaf patch groups.
upper = [g.top_core,g.top_collar,g.upper_wall];
lower = [g.bottom_core,g.bottom_collar, ...
    reshape(flipud(g.lower_wall_by_profile),1,[])];
mode = [ones(1,numel(g.top_core)+numel(g.top_collar)), ...
    2*ones(1,numel(g.upper_wall))];
permutation = zeros(1,S.npts);
upper_nodes = zeros(1,S.npts/2);
cache = struct(); count = 0;
for j = 1:numel(upper)
    a = upper(j); b = lower(j);
    ia = S.ixyzs(a):S.ixyzs(a+1)-1;
    ib = S.ixyzs(b):S.ixyzs(b+1)-1;
    key = sprintf('k%d_%d_%d',S.norders(a),S.iptype(a),mode(j));
    if ~isfield(cache,key)
        uv = S.uvs_targ(:,ia);
        if mode(j)==1
            reflected = uv([2,1],:);
        else
            reflected = [-uv(1,:);uv(2,:)];
        end
        distance = (uv(1,:).'-reflected(1,:)).^2+ ...
            (uv(2,:).'-reflected(2,:)).^2;
        [~,p] = min(distance,[],1);
        cache.(key) = p;
    end
    ib = ib(cache.(key));
    permutation(ia) = ib;
    permutation(ib) = ia;
    upper_nodes(count+(1:numel(ia))) = ia;
    count = count+numel(ia);
end
scale = max(1,max(abs(S.r),[],'all'));
reflection = [1;1;-1];
weights = S.wts(:).';
assert(count==S.npts/2 && all(permutation>0) && ...
    isequal(permutation(permutation),1:S.npts) && ...
    max(abs(S.r(:,permutation)-reflection.*S.r),[],'all')<2e4*eps(scale) && ...
    max(abs(S.n(:,permutation)-reflection.*S.n),[],'all')<2e4*eps && ...
    max(abs(weights(permutation)-weights))<2e4*eps(max(abs(weights))), ...
    'Midplane patch pairing failed.')
end
