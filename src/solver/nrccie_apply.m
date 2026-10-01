function [value,A] = nrccie_apply(S,density,T,C,opts)
%NRCCIE_APPLY Apply the NRCCIE block from surfer S to surfer T.
%   nrccie_apply(S,density,C,opts)       self, includes opts.jump * I
%   nrccie_apply(S,density,T,C,opts)     cross; T=[] also means self
% density: [ju;jv;rho], 3-by-S.npts or density(:); output follows input shape.
% C: {Cslp,Cx,Cy,Cz}, or [] for no correction. FMM also accepts a callback
% [p,gx,gy,gz]=C([Jx;Jy;Jz;rho]) for symmetry-reused corrections.
% opts: zk, alpha, fmm, eps_fmm (FMM only). All are supplied by the caller.
% opts.jump: +0.5 or -0.5, required for self; ignored for cross.
% fmm=true returns A=[]; fmm=false builds A and returns value=A*density(:).

if nargin==4
    opts = C;
    C = T;
    T = [];
end
self = isempty(T);
if self
    T = S;
end
ns = S.npts;
nt = T.npts;
k = opts.zk;
alpha = opts.alpha;
vector_input = isvector(density);
d = reshape(density,3,ns);

if opts.fmm
    % Four scalar Helmholtz potentials and their target derivatives.
    cartesian = [S.dru.*d(1,:)+S.drv.*d(2,:);d(3,:)];
    src.sources = S.r;
    src.nd = 4;
    src.charges = cartesian.*reshape(S.wts,1,ns);
    if self
        out = hfmm3d(opts.eps_fmm,k,src,2);
        p = reshape(out.pot,4,nt);
        g = reshape(out.grad,4,3,nt);
    else
        out = hfmm3d(opts.eps_fmm,k,src,0,T.r,2);
        p = reshape(out.pottarg,4,nt);
        g = reshape(out.gradtarg,4,3,nt);
    end
    gx = reshape(g(:,1,:),4,nt);
    gy = reshape(g(:,2,:),4,nt);
    gz = reshape(g(:,3,:),4,nt);
    if ~isempty(C)
        if isa(C,'function_handle')
            [cp,cx,cy,cz] = C(cartesian);
        else
            cp = (C{1}*cartesian.').';
            cx = (C{2}*cartesian.').';
            cy = (C{3}*cartesian.').';
            cz = (C{4}*cartesian.').';
        end
        p = p+cp;
        gx = gx+cx;
        gy = gy+cy;
        gz = gz+cz;
    end

    % E = ik S[J] - grad S[rho], H = curl S[J].
    E = 1i*k*p(1:3,:)-[gx(4,:);gy(4,:);gz(4,:)];
    H = [gy(3,:)-gz(2,:);gz(1,:)-gx(3,:);gx(2,:)-gy(1,:)];
    normal_e = sum(T.n.*E,1);
    tangent = -cross(T.n,H,1)+alpha*(T.n.*normal_e-E);
    value = [sum(T.dru.*tangent,1);sum(T.drv.*tangent,1); ...
        -normal_e+alpha*(gx(1,:)+gy(2,:)+gz(3,:)-1i*k*p(4,:))];
    if self
        value = value+opts.jump*d;
    end
    A = [];
else
    target = T;
    if self, target = []; end
    A = nrccie_system_block(S,target,C,opts,1:3*nt,1:3*ns);
    value = reshape(A*d(:),3,nt);
end
if vector_input
    value = value(:);
end
end
