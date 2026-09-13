function B = nrccie_quad_corr_block(S,eps_quad,zk,T,batch_size,target_ids)
% Limit near-quadrature workspace by batching target rows. Empty T means self.
if nargin<5 || isempty(batch_size), batch_size=2000; end
validateattributes(batch_size,{'numeric'},{'scalar','integer','positive'})
self=isempty(T);
if self, T=S; end
if nargin<6 || isempty(target_ids), target_ids=1:T.npts; end
nblocks=ceil(numel(target_ids)/batch_size);
rows=cell(nblocks,4);
for j=1:nblocks
    ids=target_ids((j-1)*batch_size+1:min(j*batch_size,numel(target_ids)));
    target=struct('r',T.r(:,ids));
    if self
        target.patch_id=S.patch_id(ids);
        target.uvs_targ=S.uvs_targ(:,ids);
    end
    block=cell(1,4);
    block{1}=em3d.slp.get_quad_corr_mat(S,eps_quad,zk,target);
    [block{2:4}]=em3d.sgrad.get_quad_corr_mat(S,eps_quad,zk,target);
    % Transpose chunks so sparse column pointers scale with target count.
    for k=1:4, rows{j,k}=block{k}.'; end
    clear block
end
B=cell(1,4);
for k=1:4
    B{k}=horzcat(rows{:,k}).';
    rows(:,k)={[]};
end
end
