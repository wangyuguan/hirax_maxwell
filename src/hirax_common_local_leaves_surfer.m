function [S,parts] = hirax_common_local_leaves_surfer(thickness,opts)
% Common body and four distinct tips, in TL/TR/BL/BR order; lengths in mm.
% parts.common_patch_ids / local_patch_ids index the expanded surface S.
% Default tips approximate the screenshots; opts.tip_polygons overrides them.
% Only the common portions have exact C4 symmetry.

if nargin<2 || isempty(opts), opts=struct(); end
opts=defaults(opts);
validateattributes(thickness,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.norder,{'numeric'},{'scalar','integer','>=',4});
for name={'common_rim_width','common_collar_width','cap_mesh_spacing', ...
        'local_cap_mesh_spacing','interface_spacing','surface_gap', ...
        'outline_geometry_tolerance'}
    validateattributes(opts.(name{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(opts.wall_profile_refinement,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.width_transition_y,{'numeric'},{'vector','real','finite','numel',2});
assert(numel(opts.width_transition_y)==2 && ...
    opts.width_transition_y(2)>opts.width_transition_y(1), ...
    'width_transition_y must be [lower upper] in the centered leaf frame.');

% Obtain exactly the same master outline as the original plate generator.
% Its temporary coarse cap is not used in the new geometry.
base_opts=struct('norder',4,'cap_mesh_spacing',20);
for name={'chunkie_order','chunkie_n0','chunkie_nchs','chunkie_newton_iterations'}
    if isfield(opts,name{1}), base_opts.(name{1})=opts.(name{1}); end
end
[~,base_parts]=hirax_chunkie_leaf_plate_surfer(thickness,base_opts);
[outlines,outline_info]=hirax_local_tip_outlines(base_parts,opts);
nc=numel(outline_info.common_outline_panels);

% Chunkie's narrow Gaussian rounder can have a much smaller curvature
% radius than its support width. Choose offsets from the ACTUAL curvature.
min_radius=inf;
for leaf=1:4
    for panel=1:outlines{leaf}.number_of_panels
        [~,d,~,~,curv]=curve(outlines{leaf},panel,linspace(0,1,49));
        assert(all(vecnorm(d)>1e-12),'Degenerate outline panel.');
        positive=curv>1e-10;
        if any(positive), min_radius=min(min_radius,min(1./curv(positive))); end
    end
end
if isempty(opts.tip_rim_width), opts.tip_rim_width=min(.04,.15*min_radius); end
if isempty(opts.tip_collar_width), opts.tip_collar_width=min(.04,.15*min_radius); end
validateattributes(opts.tip_rim_width,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(opts.tip_collar_width,{'numeric'},{'scalar','real','finite','positive'});

% Refinement depends only on master geometry, never on surface norder.
% Each common panel is refined identically in all four outlines.
for leaf=1:4
    [outlines{leaf},counts]=refine_outline(outlines{leaf},opts);
    ncommon(leaf)=sum(counts(1:nc)); %#ok<AGROW>
end
assert(all(ncommon==ncommon(1)),'Common panel refinement must agree.');
nc=ncommon(1);
for leaf=2:4
    assert(isequal(outlines{1}.position_coefficients(:,:,1:nc), ...
        outlines{leaf}.position_coefficients(:,:,1:nc)), ...
        'Common outline coefficients must be exactly identical.');
end

% Common core: one polygon along the common arc, closed by a fixed chord.
% Local core: the remaining boundary plus the reversed SAME chord nodes.
core_vertices=cell(1,4);
for leaf=1:4
    ol=outlines{leaf}; v=zeros(2,ol.number_of_panels);
    for panel=1:ol.number_of_panels
        [~,~,v(:,panel)]=offset_curve(ol,panel,0,opts);
    end
    core_vertices{leaf}=v;
end
[~,~,left_attach]=offset_curve(outlines{1},nc,1,opts);
right_attach=core_vertices{1}(:,1);
for leaf=1:4
    % Pin both sides before constructing either piece of the cap mesh.
    core_vertices{leaf}(:,1)=right_attach;
    core_vertices{leaf}(:,nc+1)=left_attach;
end
nt=ceil(norm(left_attach-right_attach)/opts.interface_spacing);
interface=left_attach+(right_attach-left_attach).*linspace(0,1,nt+1);
common_polygon=[core_vertices{1}(:,1:nc),interface(:,1:end-1)];
[common_nodes,common_faces]=triangulate_polygon(common_polygon,opts.cap_mesh_spacing);
[common,common_groups]=surface_piece(outlines{1},1:nc,core_vertices{1}, ...
    common_nodes,common_faces,thickness,opts);
unique_surfaces=cell(1,4); local_core=cell(1,4); local_groups=cell(1,4);
for leaf=1:4
    local_polygon=[core_vertices{leaf}(:,nc+1:end),right_attach, ...
        interface(:,end-1:-1:2)];
    [nodes,faces]=triangulate_polygon(local_polygon,opts.local_cap_mesh_spacing);
    [unique_surfaces{leaf},local_groups{leaf}]=surface_piece(outlines{leaf}, ...
        nc+1:outlines{leaf}.number_of_panels,core_vertices{leaf}, ...
        nodes,faces,thickness,opts);
    local_core{leaf}=struct('nodes',nodes,'faces',faces,'polygon',local_polygon);
end

parts=struct();
parts.local_common_surfer=common;
parts.local_unique_surfers=unique_surfaces;
parts.common_core=struct('nodes',common_nodes,'faces',common_faces, ...
    'polygon',common_polygon);
parts.local_core=local_core;
parts.interface_nodes=interface;
parts.common_groups=common_groups;
parts.local_groups=local_groups;
parts.outlines=outlines;
parts.outline_info=outline_info;
parts.outline_info.indexing='Outline-helper indices before cap/rim refinement';
parts.common_outline_panels=1:nc;
parts.local_outline_panels=cellfun(@(o) nc+1:o.number_of_panels, ...
    outlines,'UniformOutput',false);
parts.thickness=thickness;
parts.options=opts;
parts.minimum_outline_convex_radius=min_radius;
parts.is_full_fourfold_symmetric=false;
parts.leaf_names={'upper_left','upper_right','lower_left','lower_right'};
parts.design_note=['Image-informed geometric prototype; central dimensions ' ...
    'are editable estimates, not a Gerber/CAD reconstruction.'];
if ~outline_info.approximate_image_informed_defaults
    parts.design_note='Shared original master body with caller-supplied central waypoint chains.';
end
parts.leaves=cell(1,4);
parts.common_patch_ids=cell(1,4); parts.local_patch_ids=cell(1,4);
parts.common_node_ids=cell(1,4); parts.local_node_ids=cell(1,4);
angles=[pi/4,-pi/4,3*pi/4,-3*pi/4];
patch_offset=0; node_offset=0;
% Anchor the unchanged long flanks, NOT the tips' extrema. The raw master
% flanks satisfy v=abs(u)+2.83, corresponding to a 2.83*sqrt(2) mm gap.
center=base_parts.design.leaf_center_raw;
center(2)=center(2)+(opts.surface_gap-2.83*sqrt(2))/sqrt(2);
for leaf=1:4
    a=angles(leaf); R=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    shift=R*[center;0];
    local_leaf=merge([common,unique_surfaces{leaf}]);
    parts.leaves{leaf}=affine_transf(local_leaf,R,shift);
    parts.rotation(:,:,leaf)=R;
    parts.shift(:,leaf)=shift;
    parts.common_patch_ids{leaf}=patch_offset+(1:common.npatches);
    parts.local_patch_ids{leaf}=patch_offset+common.npatches+ ...
        (1:unique_surfaces{leaf}.npatches);
    parts.common_node_ids{leaf}=node_offset+(1:common.npts);
    parts.local_node_ids{leaf}=node_offset+common.npts+(1:unique_surfaces{leaf}.npts);
    patch_offset=patch_offset+local_leaf.npatches;
    node_offset=node_offset+local_leaf.npts;
end
S=merge([parts.leaves{:}]);
end

function o=defaults(o)
d=struct('norder',6,'common_rim_width',.15,'common_collar_width',.6, ...
    'tip_rim_width',[],'tip_collar_width',[],'width_transition_y',[-55,-40], ...
    'cap_mesh_spacing',6,'local_cap_mesh_spacing',2,'interface_spacing',2, ...
    'wall_profile_refinement',3,'surface_gap',4,'outline_geometry_tolerance',1e-6);
names=fieldnames(d);
for j=1:numel(names)
    if ~isfield(o,names{j}), o.(names{j})=d.(names{j}); end
end
end

function [r,d,n,dn,kappa]=curve(o,p,q)
% Legendre derivative coefficients, with q in [0,1].
k=o.order; D=zeros(k);
for col=2:k
    rows=col-1:-2:1; D(rows,col)=2*rows-1;
end
P=reshape(lege.pols(2*q-1,k-1),k,[]);
c=o.position_coefficients(:,:,p);
r=c*P; d=2*(c*D.')*P; dd=4*(c*D.'*D.')*P;
speed=vecnorm(d); tang=d./speed;
n=[tang(2,:);-tang(1,:)];
dt=(dd-tang.*sum(tang.*dd,1))./speed;
dn=[dt(2,:);-dt(1,:)];
kappa=(d(1,:).*dd(2,:)-d(2,:).*dd(1,:))./speed.^3;
end

function [w,dw]=width(r,d,lo,hi,opts)
a=opts.width_transition_y(1); b=opts.width_transition_y(2);
t=max(0,min(1,(r(2,:)-a)/(b-a)));
s=t.^3.*(10-15*t+6*t.^2);
ds=30*t.^2.*(1-t).^2/(b-a);
w=lo+(hi-lo)*s;
dw=(hi-lo)*ds.*d(2,:);
end

function [inner,di,core,dc,r,d,n]=offset_curve(o,p,q,opts)
[r,d,n,dn,kappa]=curve(o,p,q);
[rw,drw]=width(r,d,opts.tip_rim_width,opts.common_rim_width,opts);
[cw,dcw]=width(r,d,opts.tip_collar_width,opts.common_collar_width,opts);
assert(all(1-(rw+cw).*kappa>.1), ...
    'hirax:OffsetFold','Rim/collar too wide for a corner; reduce widths or increase rounding.');
inner=r-rw.*n; di=d-drw.*n-rw.*dn;
core=r-(rw+cw).*n; dc=d-(drw+dcw).*n-(rw+cw).*dn;
end

function [out,counts]=refine_outline(o,opts)
coeff=cell(1,o.number_of_panels); counts=zeros(1,o.number_of_panels);
[x,~,V]=lege.exps(o.order); q=(x(:).'+1)/2;
[xcheck,~,Vcheck]=lege.exps(5);
qcheck=(xcheck(:).'+1)/2;
Pcheck=reshape(lege.pols(linspace(-1,1,33),4),5,[]);
for p=1:o.number_of_panels
    intervals=[0,1]; accepted=zeros(0,2);
    while ~isempty(intervals)
        ab=intervals(1,:); intervals(1,:)=[];
        tt=linspace(ab(1),ab(2),33);
        [inner,di,cc,~,rr,~,normals]=offset_curve(o,p,tt,opts);
        ends=cc(:,[1,end]); s=linspace(0,1,33);
        chord=ends(:,1)+(ends(:,2)-ends(:,1)).*s;
        dev=max(vecnorm(cc-chord));
        minwidth=min(vecnorm(cc-inner));
        tangchange=max(vecnorm(normals-normals(:,1)));
        a=di*(ab(2)-ab(1)); b=chord-inner;
        c=repmat(ends(:,2)-ends(:,1),1,numel(s));
        jac0=a(1,:).*b(2,:)-a(2,:).*b(1,:);
        jac1=c(1,:).*b(2,:)-c(2,:).*b(1,:);
        % Resolve even the lowest supported order on the SAME panelization.
        % Normal offsets can vary faster than the original position curve
        % at a tightly rounded corner; checking sagitta alone misses this.
        qt=ab(1)+(ab(2)-ab(1))*qcheck;
        [ii,~,ci,~,ri]=offset_curve(o,p,qt,opts);
        values=[ri;ii;ci]; approx=(Vcheck*values.').'*Pcheck;
        approximation_error=max(abs(approx-[rr;inner;cc]),[],'all');
        good=dev<.2*minwidth && tangchange<.3 && all(jac0>0) && all(jac1>0) ...
            && approximation_error<opts.outline_geometry_tolerance;
        if good
            accepted(end+1,:)=ab; %#ok<AGROW>
        else
            assert(ab(2)-ab(1)>2^-16,'Could not produce an unfolded cap collar.');
            mid=mean(ab); intervals=[ab(1),mid;mid,ab(2);intervals]; %#ok<AGROW>
        end
    end
    accepted=sortrows(accepted,1); counts(p)=size(accepted,1);
    cp=zeros(2,o.order,counts(p));
    for j=1:counts(p)
        ab=accepted(j,:);
        if all(ab==[0,1])
            cp(:,:,j)=o.position_coefficients(:,:,p);
        else
            rr=curve(o,p,ab(1)+(ab(2)-ab(1))*q);
            cp(:,:,j)=(V*rr.').';
        end
    end
    coeff{p}=cp;
end
out=o; out.position_coefficients=cat(3,coeff{:});
out.number_of_panels=size(out.position_coefficients,3);
out.panel_starts=zeros(2,out.number_of_panels);
out.panel_start_derivatives=out.panel_starts;
for p=1:out.number_of_panels
    [out.panel_starts(:,p),out.panel_start_derivatives(:,p)]=curve(out,p,0);
end
end

function [nodes,faces]=triangulate_polygon(vertices,spacing)
% Constraints explicitly include every shared interface node; never refine
% this boundary independently in the two cap subdomains.
delta=vertices(:,[2:end,1])-vertices;
assert(all(vecnorm(delta)>1e-10),'Repeated polygon boundary vertices.');
[xx,yy]=meshgrid(min(vertices(1,:)):spacing:max(vertices(1,:)), ...
    min(vertices(2,:)):spacing:max(vertices(2,:)));
keep=inpolygon(xx(:),yy(:),vertices(1,:),vertices(2,:));
interior=[xx(keep).';yy(keep).'];
for edge=1:size(vertices,2)
    a=vertices(:,edge); v=delta(:,edge);
    t=max(0,min(1,sum((interior-a).*v,1)/sum(v.^2)));
    interior=interior(:,vecnorm(interior-(a+v.*t))>.2*spacing);
end
points=[vertices,interior].'; nv=size(vertices,2);
dt=delaunayTriangulation(points,[(1:nv).',[2:nv,1].']);
nodes=dt.Points.'; faces=dt.ConnectivityList(isInterior(dt),:).';
a=nodes(:,faces(2,:))-nodes(:,faces(1,:));
b=nodes(:,faces(3,:))-nodes(:,faces(1,:));
area=a(1,:).*b(2,:)-a(2,:).*b(1,:);
assert(all(abs(area)>1e-14),'Degenerate cap triangle.');
faces([2,3],area<0)=faces([3,2],area<0);
end

function [S,g]=surface_piece(o,panels,vertices,nodes,faces,thickness,opts)
p=opts.norder; uv=koorn.rv_nodes(p); nt=size(uv,2);
q=(polytens.lege.nodes(p)+1)/2; nq=size(q,2);
nf=size(faces,2); np=numel(panels); nw=opts.wall_profile_refinement;
top=zeros(12,nt*nf); bottom=top;
for j=1:nf
    ids=(j-1)*nt+(1:nt); vv=nodes(:,faces(:,j));
    top(:,ids)=triangle(vv,thickness/2,uv);
    bottom(:,ids)=triangle(vv(:,[1,3,2]),-thickness/2,uv);
end
tc=zeros(12,nq*np); bc=tc;
uw=zeros(12,nq*np*nw); lw=uw;
for j=1:np
    panel=panels(j); ids=(j-1)*nq+(1:nq);
    ends=vertices(:,[panel,mod(panel,o.number_of_panels)+1]);
    tc(:,ids)=collar(o,panel,ends,q,thickness/2,1,opts);
    bc(:,ids)=collar(o,panel,ends,q,-thickness/2,-1,opts);
    for h=1:nw
        ids=((j-1)*nw+h-1)*nq+(1:nq); interval=(h-1:h)/nw;
        uw(:,ids)=wall(o,panel,q,interval,thickness/2,1,opts);
        lw(:,ids)=wall(o,panel,q,interval,thickness/2,-1,opts);
    end
end
objects=[surfer(nf,p,top,1),surfer(np,p,tc,11), ...
    surfer(np*nw,p,uw,11),surfer(np*nw,p,lw,11), ...
    surfer(np,p,bc,11),surfer(nf,p,bottom,1)];
S=merge(objects);
names={'top_core','top_collar','upper_wall','lower_wall','bottom_collar','bottom_core'};
offset=0; g=struct();
for j=1:numel(objects)
    g.(names{j})=offset+(1:objects(j).npatches); offset=offset+objects(j).npatches;
end
g.upper_wall_by_profile=reshape(g.upper_wall,nw,np);
g.lower_wall_by_profile=reshape(g.lower_wall,nw,np);
g.outline_panels=panels;
end

function v=triangle(x,z,uv)
n=size(uv,2); a=x(:,2)-x(:,1); b=x(:,3)-x(:,1);
r=[x(:,1)+a.*uv(1,:)+b.*uv(2,:);z*ones(1,n)];
v=pack(r,repmat([a;0],1,n),repmat([b;0],1,n));
end

function v=collar(o,p,ends,q,z,orientation,opts)
if orientation>0, t=q(1,:); s=q(2,:); else, t=q(2,:); s=q(1,:); end
[inner,di]=offset_curve(o,p,t,opts);
line=ends(:,1)+(ends(:,2)-ends(:,1)).*t;
dt=(1-s).*di+s.*(ends(:,2)-ends(:,1)); ds=line-inner;
assert(all(dt(1,:).*ds(2,:)-dt(2,:).*ds(1,:)>0),'Folded cap collar.');
r=[(1-s).*inner+s.*line;z*ones(size(t))];
dt=.5*[dt;zeros(size(t))]; ds=.5*[ds;zeros(size(t))];
if orientation>0, v=pack(r,dt,ds); else, v=pack(r,ds,dt); end
end

function v=wall(o,p,q,interval,halfthickness,side,opts)
[inner,di,~,~,outer,doo]=offset_curve(o,p,q(2,:),opts);
span=diff(interval); s=interval(1)+span*q(1,:);
if side<0, s=1-s; span=-span; end
[rad,dr]=bezier([0,.4,.8,1,1,1],s);
[height,dh]=bezier([1,1,1,.8,.4,0],s);
r=[inner+rad.*(outer-inner);side*halfthickness*height];
du=.5*span*[(outer-inner).*dr;side*halfthickness*dh];
dv=.5*[di+rad.*(doo-di);zeros(size(s))];
v=pack(r,du,dv);
end

function [value,derivative]=bezier(control,s)
degree=numel(control)-1; value=zeros(size(s)); derivative=value;
for j=0:degree
    value=value+control(j+1)*nchoosek(degree,j)*(1-s).^(degree-j).*s.^j;
end
dc=degree*diff(control);
for j=0:degree-1
    derivative=derivative+dc(j+1)*nchoosek(degree-1,j)*(1-s).^(degree-1-j).*s.^j;
end
end

function v=pack(r,du,dv)
n=cross(du,dv,1); jac=vecnorm(n);
assert(all(jac>0) && all(isfinite(jac)),'Degenerate surface patch.');
v=[r;du;dv;n./jac];
end
