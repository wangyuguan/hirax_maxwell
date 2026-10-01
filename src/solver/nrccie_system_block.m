function A = nrccie_system_block(S,T,C,opts,rows,cols)
%NRCCIE_SYSTEM_BLOCK Selected entries of the corrected NRCCIE system.
% rows/cols are global interleaved [ju;jv;rho] DOF indices on T/S.
% T=[] means self and includes opts.jump on matching global DOFs.
% Explicit T means cross, with no jump. C contains the four scalar
% correction matrices for the full target/source surfaces, or is empty.
self = isempty(T);
if self, T = S; end
rows = rows(:).';
cols = cols(:).';
A = complex(zeros(numel(rows),numel(cols)));
if isempty(rows) || isempty(cols), return; end

[target_nodes,~,target_map] = unique(ceil(rows/3));
[source_nodes,~,source_map] = unique(ceil(cols/3));
target_map = target_map(:);
ns = numel(source_nodes);
k = opts.zk;
alpha = opts.alpha;
weights = reshape(S.wts(source_nodes),1,ns);

% A current column carries its source tangent; a charge column carries zero.
% Keep only the requested scalar DOFs, including their ordering/repetitions.
source_current = zeros(3,numel(cols));
component = mod(cols-1,3)+1;
source_current(:,component==1) = S.dru(:,ceil(cols(component==1)/3));
source_current(:,component==2) = S.drv(:,ceil(cols(component==2)/3));
charge = component==3;
sx = source_current(1,:);
sy = source_current(2,:);
sz = source_current(3,:);

% Bound kernel workspace even when FLAM requests a long off-diagonal block.
% No full 3*nt-by-3*ns block or Cartesian E/H matrices are formed.
for first = 1:32:numel(target_nodes)
    last = min(first+31,numel(target_nodes));
    nodes = target_nodes(first:last);
    dx = T.r(1,nodes).'-S.r(1,source_nodes);
    dy = T.r(2,nodes).'-S.r(2,source_nodes);
    dz = T.r(3,nodes).'-S.r(3,source_nodes);
    r = sqrt(dx.^2+dy.^2+dz.^2);
    rinv = 1./r;
    if self, rinv(nodes(:)==source_nodes) = 0; end
    green = exp(1i*k*r).*rinv.*weights/(4*pi);
    gradient = green.*(1i*k*r-1).*rinv.^2;
    gx = dx.*gradient;
    gy = dy.*gradient;
    gz = dz.*gradient;
    if ~isempty(C)
        green = green+C{1}(nodes,source_nodes);
        gx = gx+C{2}(nodes,source_nodes);
        gy = gy+C{3}(nodes,source_nodes);
        gz = gz+C{4}(nodes,source_nodes);
    end
    green = green(:,source_map);
    gx = gx(:,source_map);
    gy = gy(:,source_map);
    gz = gz(:,source_map);
    n = T.n(:,nodes).';
    normal_g = n(:,1).*gx+n(:,2).*gy+n(:,3).*gz;
    normal_current = n(:,1).*sx+n(:,2).*sy+n(:,3).*sz;

    for c = 1:3
        selected = find(target_map>=first & target_map<=last ...
            & mod(rows(:)-1,3)+1==c);
        if isempty(selected), continue; end
        if c==3
            values = -1i*k*green.*normal_current ...
                +alpha*(gx.*sx+gy.*sy+gz.*sz);
            values(:,charge) = normal_g(:,charge)-alpha*1i*k*green(:,charge);
        else
            if c==1, tangent = T.dru(:,nodes).';
            else, tangent = T.drv(:,nodes).'; end
            tangent_g = tangent(:,1).*gx+tangent(:,2).*gy+tangent(:,3).*gz;
            tangent_current = tangent(:,1).*sx+tangent(:,2).*sy+tangent(:,3).*sz;
            tangent_normal = sum(tangent.*n,2);
            values = tangent_current.*normal_g-tangent_g.*normal_current ...
                +alpha*1i*k*green.*(tangent_normal.*normal_current-tangent_current);
            values(:,charge) = alpha*(tangent_g(:,charge) ...
                -tangent_normal.*normal_g(:,charge));
        end
        if self, values = values+opts.jump*(3*nodes(:)-3+c==cols); end
        A(selected,:) = values(target_map(selected)-first+1,:);
    end
end
end
