function value = nrccie_apply(S,density,C,opts)
%
%  nrccie_apply
%    Apply the corrected NRCCIE operator using the FMM.
%
%  Input arguments:
%    * S: source and target surfer
%    * density: interleaved [ju,jv,rho], as a vector, a (3,npts) array,
%        or a (3*npts,nrhs) array for several right-hand sides
%    * C: quadrature corrections {Cslp,Cx,Cy,Cz}, [] for no correction,
%        or callback [p,gx,gy,gz] = C([Jx;Jy;Jz;rho])
%    * opts: zk, alpha, jump, and eps_fmm
%
%  Output arguments:
%    * value: column vector, (3,npts) array, or (3*npts,nrhs) array,
%        matching the input density representation
%
    ns = S.npts;
    nt = ns;
    k = opts.zk;
    alpha = opts.alpha;
    vector_input = isvector(density);
    surface_input = ~vector_input && isequal(size(density),[3,ns]);
    if vector_input || surface_input
      assert(numel(density)==3*ns,'Density size does not match the source.')
      density_columns = reshape(density,3*ns,1);
    else
      assert(ismatrix(density) && size(density,1)==3*ns, ...
          'A density matrix must have 3*S.npts rows.')
      density_columns = density;
    end
    nrhs = size(density_columns,2);
    d = permute(reshape(density_columns,3,ns,nrhs),[1,3,2]);

    % Four scalar Helmholtz potentials and their target derivatives.
    cartesian = zeros(4,nrhs,ns,'like',density_columns);
    cartesian(1:3,:,:) = reshape(S.dru,3,1,ns).*d(1,:,:) + ...
        reshape(S.drv,3,1,ns).*d(2,:,:);
    cartesian(4,:,:) = d(3,:,:);
    cartesian = reshape(cartesian,4*nrhs,ns);
    src.sources = S.r;
    src.nd = 4*nrhs;
    src.charges = cartesian.*reshape(S.wts,1,ns);
    out = hfmm3d(opts.eps_fmm,k,src,2);
    p = reshape(out.pot,4*nrhs,nt);
    g = reshape(out.grad,4*nrhs,3,nt);
    gx = reshape(g(:,1,:),4*nrhs,nt);
    gy = reshape(g(:,2,:),4*nrhs,nt);
    gz = reshape(g(:,3,:),4*nrhs,nt);
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

    p = reshape(p,4,nrhs,nt);
    gx = reshape(gx,4,nrhs,nt);
    gy = reshape(gy,4,nrhs,nt);
    gz = reshape(gz,4,nrhs,nt);

    % E = ik S[J] - grad S[rho], H = curl S[J].
    E = 1i*k*p(1:3,:,:)-cat(1,gx(4,:,:),gy(4,:,:),gz(4,:,:));
    H = cat(1,gy(3,:,:)-gz(2,:,:),gz(1,:,:)-gx(3,:,:), ...
        gx(2,:,:)-gy(1,:,:));
    normal = repmat(reshape(S.n,3,1,nt),1,nrhs,1);
    normal_e = sum(normal.*E,1);
    tangent = -cross(normal,H,1)+alpha*(normal.*normal_e-E);
    value = cat(1,sum(reshape(S.dru,3,1,nt).*tangent,1), ...
        sum(reshape(S.drv,3,1,nt).*tangent,1), ...
        -normal_e+alpha*(gx(1,:,:)+gy(2,:,:)+gz(3,:,:)-1i*k*p(4,:,:)));
    value = value+opts.jump*d;
    columns = reshape(permute(value,[1,3,2]),3*nt,nrhs);
    if vector_input
      value = columns(:,1);
    elseif surface_input
      value = reshape(columns(:,1),3,nt);
    else
      value = columns;
    end
end
