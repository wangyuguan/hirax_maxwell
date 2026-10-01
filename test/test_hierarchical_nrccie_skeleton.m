% Small check of both hierarchical maps; no integral-equation solve or plots.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(genpath(fullfile(root,'src')))
clear pth dir
rng(37)

% Chebyshev patches include their corners, so shared edges match exactly.
S = geometries.ellipsoid([1.2,0.85,0.65],[1,1,1],[0;0;0],5,12);
T = geometries.ellipsoid(4*[1,1,1],[1,1,1],[0;0;0],6,11);
zk = 0.9;
alpha = 1;
tol = 1e-8;
max_patches = 2;
rotations = sample_so3(128);
validation_rotations = sample_so3(32);
assert(min(vecnorm(T.r))>max(vecnorm(S.r)))
source_w = sqrt(repelem(S.wts(:),3));
proxy_w = sqrt(repelem(T.wts(:),6));

[Jin,blocks_in,Pmerge_in,info_in] = hierarchical_nrccie_incoming_skeleton( ...
    S,zk,alpha,rotations,tol,max_patches);
[Jout,blocks_out,Pmerge_out,info_out] = hierarchical_nrccie_outgoing_skeleton( ...
    S,T,zk,tol,max_patches);

% Every node occurs once, and every part connects through complete edges.
assert(isequal(info_in.part_id,info_out.part_id))
for incoming = [true,false]
    if incoming
        blocks = blocks_in;
        info = info_in;
        J = Jin;
    else
        blocks = blocks_out;
        info = info_out;
        J = Jout;
    end
    assert(numel(blocks)>1)
    assert(all(info.part_patch_counts<=max_patches))
    assert(isequal(sort([blocks.nodes]),1:S.npts))
    assert(isequal([blocks.skeleton_nodes],info.candidate_nodes))
    assert(isequal([blocks.candidate_dofs],1:3*numel(info.candidate_nodes)))
    skeleton_nodes = (J(1:3:end)+2)/3;
    assert(isequal(reshape(J,3,[]),3*skeleton_nodes+(-2:0).'))
    assert(numel(unique(skeleton_nodes))==numel(skeleton_nodes))
    for part = 1:numel(blocks)
        patches = find(info.part_id==part);
        assert(all(conncomp(graph(info.adjacency(patches,patches)))==1))
        assert(all(ismember(S.patch_id(blocks(part).nodes),patches)))
        assert(all(ismember(blocks(part).skeleton_nodes,blocks(part).nodes)))
    end
end

% Full references are deliberately small and are used only by this test.
B = nrccie_planewave_rhs(S,zk,alpha,rotations);
Btest = nrccie_planewave_rhs(S,zk,alpha,validation_rotations);
Pin = complex(zeros(3*S.npts,numel(Jin)));
factored_rhs = complex(zeros(size(Btest)));
candidate_rhs = Pmerge_in*Btest(Jin,:);
for part = 1:numel(blocks_in)
    block = blocks_in(part);
    dofs = reshape(3*block.nodes+(-2:0).',1,[]);
    Pin(dofs,:) = block.P*Pmerge_in(block.candidate_dofs,:);
    factored_rhs(dofs,:) = block.P*candidate_rhs(block.candidate_dofs,:);
end
assert(norm(Pin(Jin,:)-eye(numel(Jin)),'fro')<1e-11)
assert(norm(factored_rhs-Pin*Btest(Jin,:),'fro')<1e-12*norm(Btest,'fro'))
training_error = norm(source_w.*(B-Pin*B(Jin,:)),'fro') ...
    /norm(source_w.*B,'fro');
heldout_residual = source_w.*(Btest-factored_rhs);
incoming_error = norm(heldout_residual,'fro')/norm(source_w.*Btest,'fro');
worst_wave_error = max(vecnorm(heldout_residual)./vecnorm(source_w.*Btest));
assert(training_error<1.1*tol)
assert(training_error<=1.1*info_in.training_error_bound+200*eps)
assert(incoming_error<3*tol && worst_wave_error<10*tol)

K = nrccie_matrix(S,T,zk);
target_rotations = sample_so3(64);
validation_targets.r = 8*reshape(target_rotations(:,3,:),3,[]);
Kvalidation = nrccie_matrix(S,validation_targets,zk);
Pout = complex(zeros(numel(Jout),3*S.npts));
q = randn(3*S.npts,4)+1i*randn(3*S.npts,4);
density = (q./vecnorm(q))./source_w;
equivalent_density = complex(zeros(numel(Jout),size(density,2)));
for part = 1:numel(blocks_out)
    block = blocks_out(part);
    dofs = reshape(3*block.nodes+(-2:0).',1,[]);
    merge_part = Pmerge_out(:,block.candidate_dofs);
    Pout(:,dofs) = merge_part*block.P;
    equivalent_density = equivalent_density+merge_part*(block.P*density(dofs,:));
end
assert(norm(Pout(:,Jout)-eye(numel(Jout)),'fro')<1e-11)
assert(norm(equivalent_density-Pout*density,'fro') ...
    <1e-12*norm(Pout*density,'fro'))
proxy_error = norm((proxy_w.*(K-K(:,Jout)*Pout))./source_w.','fro') ...
    /norm((proxy_w.*K)./source_w.','fro');
validation_error = norm((Kvalidation-Kvalidation(:,Jout)*Pout)./source_w.','fro') ...
    /norm(Kvalidation./source_w.','fro');
density_error = norm(Kvalidation*density-Kvalidation(:,Jout)*equivalent_density,'fro') ...
    /norm(Kvalidation*density,'fro');
assert(proxy_error<1.1*tol)
assert(proxy_error<=1.1*info_out.proxy_error_bound+200*eps)
assert(validation_error<3*tol && density_error<10*tol)

fprintf('Surface: %d nodes in %d connected parts; proxy: %d nodes\n', ...
    S.npts,numel(blocks_in),T.npts)
fprintf('Incoming: %d candidates -> %d complete nodes; outgoing: %d -> %d\n', ...
    numel(info_in.candidate_nodes),numel(Jin)/3, ...
    numel(info_out.candidate_nodes),numel(Jout)/3)
fprintf('Incoming errors: training %.3e (bound %.3e), held-out %.3e, worst %.3e\n', ...
    training_error,info_in.training_error_bound,incoming_error,worst_wave_error)
fprintf('Outgoing errors: proxy %.3e (bound %.3e), exterior field %.3e, density %.3e\n', ...
    proxy_error,info_out.proxy_error_bound,validation_error,density_error)
disp('HIERARCHICAL_NRCCIE_SKELETON_PASS')
