function proxyfun = nrccie_rskelf_proxy(S,C,opts,tol)
%NRCCIE_RSKELF_PROXY Local, two-sided proxy callback for FLAM rskelf.
% Proxy points lie on ctr + 1.3*l.*unit_sphere, enclosing the current box
% inside its geometric neighbor region. Those neighbors remain explicit.
% The four Helmholtz potentials SJx,SJy,SJz,Srho determine the full NRCCIE
% operator, including div SJ-ik Srho. Sample their outgoing values and the
% adjoint incoming operator with four independent Cartesian source channels.
% This uses exterior Helmholtz uniqueness (real positive zk), not a density
% restricted to two tangents on the artificial surface.
% opts.proxy_order is the polar order; there are 2*order^2 proxy points.
% Increasing it controls proxy sampling error separately from FLAM's ID tol.

order = max(12,ceil(-2*log10(tol)));
if isfield(opts,'proxy_order'), order = opts.proxy_order; end
assert(order>=2 && order==fix(order),'proxy_order must be an integer >= 2.')
% Account for waves on the largest proxy-eligible box; a conservative root
% size keeps the same reference sampling available at every tree level.
if ~isfield(opts,'proxy_order')
    order = order+ceil(abs(opts.zk)*max(max(S.r,[],2)-min(S.r,[],2)));
end
[z,w] = glegquad(order);
phi = (0:2*order-1)*pi/order;
unit = [reshape(sqrt(1-z.^2).*cos(phi),1,[]); ...
        reshape(sqrt(1-z.^2).*sin(phi),1,[]); ...
        repmat(z(:).',1,2*order)];
row_scale = sqrt(S.npts*repmat(w(:),2*order,1)/(4*order));
mean_weight = mean(S.wts);
correction_support = sparse(S.npts,S.npts);
for j = 1:numel(C)
    correction_support = correction_support+spones(C{j});
end
correction_support = spones(correction_support+correction_support.');
matrix_block = @(rows,cols)nrccie_system_block(S,[],C,opts,rows,cols);
proxyfun = @proxy_block;

    function [K,nbr] = proxy_block(~,slf,nbr,l,ctr)
        slf = slf(:).';
        points = ctr(:)+1.3*l(:).*unit;
        np = size(points,2);
        nodes = ceil(slf/3);
        component = mod(slf-1,3)+1;
        charge = component==3;
        tangent = zeros(3,numel(slf));
        tangent(:,component==1) = S.dru(:,nodes(component==1));
        tangent(:,component==2) = S.drv(:,nodes(component==2));
        normal = S.n(:,nodes);
        tangent_normal = sum(tangent.*normal,1);

        % Gradient with respect to the local target, for incoming samples.
        dx = S.r(1,nodes)-points(1,:).';
        dy = S.r(2,nodes)-points(2,:).';
        dz = S.r(3,nodes)-points(3,:).';
        r = sqrt(dx.^2+dy.^2+dz.^2);
        green = exp(1i*opts.zk*r)./(4*pi*r);
        derivative = green.*(1i*opts.zk*r-1)./r.^2;
        gradient = {dx.*derivative,dy.*derivative,dz.*derivative};
        normal_g = gradient{1}.*normal(1,:)+gradient{2}.*normal(2,:) ...
            +gradient{3}.*normal(3,:);
        tangent_g = gradient{1}.*tangent(1,:)+gradient{2}.*tangent(2,:) ...
            +gradient{3}.*tangent(3,:);
        ik_green = 1i*opts.zk*green;
        alpha = opts.alpha;
        % Put potential traces on the differentiated operator's scale.
        outgoing = green.*reshape(S.wts(nodes),1,[]) ...
            *(abs(opts.zk)+1/min(l))*max(1,abs(alpha));
        K = complex(zeros(8*np,numel(slf)));
        for c = 1:3
            rows = (c-1)*np+(1:np);
            K(rows,:) = row_scale.*outgoing.*tangent(c,:);
            incoming = tangent(c,:).*normal_g-normal(c,:).*tangent_g ...
                +alpha*ik_green.*(tangent_normal.*normal(c,:)-tangent(c,:));
            incoming(:,charge) = -ik_green(:,charge).*normal(c,charge) ...
                +alpha*gradient{c}(:,charge);
            K(4*np+rows,:) = row_scale.*conj(mean_weight*incoming);
        end
        K(3*np+(1:np),:) = row_scale.*outgoing.*charge;
        incoming = alpha*(tangent_g-tangent_normal.*normal_g);
        incoming(:,charge) = normal_g(:,charge)-alpha*ik_green(:,charge);
        K(7*np+(1:np),:) = row_scale.*conj(mean_weight*incoming);

        % Corrections need not fit inside FLAM's geometric neighbor region.
        % Sample their full interactions explicitly, in BOTH directions.
        % These are sample rows, not added nbr indices: some original DOFs
        % may already have been eliminated and must not re-enter the tree.
        related = find(any(correction_support(:,unique(nodes)),2));
        inside = all(abs((S.r(:,related)-ctr(:))./l(:))<=0.5+16*eps,1);
        related = setdiff(related(~inside),unique(ceil(nbr/3)));
        extra = setdiff(reshape(3*related(:).'+(-2:0).',1,[]),[slf(:);nbr(:)]);
        if ~isempty(extra)
            K = [K;matrix_block(extra,slf);matrix_block(slf,extra)'];
        end
    end
end
