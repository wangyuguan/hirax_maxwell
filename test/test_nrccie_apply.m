% Small regression test for modular NRCCIE source/target application.
clear
test_dir = fileparts(mfilename('fullpath'));
root = fileparts(test_dir);
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
addpath(fullfile(root,'..','fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(root,'src')))
clear pth dir

rng(19)
zk = 0.7;
alpha = 0.8;
eps_quad = 1e-10;
S = geometries.sphere(1,1,[0;0;0],2,1);
T = geometries.sphere(1.15,1,[0;0;0],3,1);
density = randn(3,S.npts)+1i*randn(3,S.npts);

fprintf('Building small self and source-to-target corrections ...\n')
Cself = nrccie_quad_corr_block(S,eps_quad,zk);
Ccross = nrccie_quad_corr_block(S,eps_quad,zk,T);
assert(sum(cellfun(@nnz,Cself))>0)
assert(sum(cellfun(@nnz,Ccross))>0)
dense_options = struct('zk',zk,'alpha',alpha,'fmm',false,'jump',0.5);
fmm_options = struct('zk',zk,'alpha',alpha, ...
    'fmm',true,'eps_fmm',1e-11,'jump',0.5);

%% Self-to-self: explicit matrix and FMM agree, and short call adds jump.
[self_dense,Aself] = nrccie_apply(S,density,Cself,dense_options);
assert(isequal(size(self_dense),[3,S.npts]))
assert(isequal(size(Aself),[3*S.npts,3*S.npts]))
assert(relative_error(self_dense(:),Aself*density(:))<5e-14)

[self_fmm,Afmm] = nrccie_apply(S,density,[],Cself,fmm_options);
assert(isempty(Afmm))
assert(relative_error(self_fmm,self_dense)<5e-8)

correction_apply = @(d)apply_full_correction(Cself,d);
self_callback = nrccie_apply(S,density,correction_apply,fmm_options);
assert(relative_error(self_callback,self_fmm)<5e-13)

%% Changing +1/2 I to -1/2 I subtracts exactly one identity.
negative_dense = dense_options;
negative_dense.jump = -0.5;
negative_fmm = fmm_options;
negative_fmm.jump = -0.5;
[self_minus,Aminus] = nrccie_apply(S,density,Cself,negative_dense);
self_minus_fmm = nrccie_apply(S,density,Cself,negative_fmm);
assert(relative_error(self_minus,self_dense-density)<5e-14)
assert(norm(Aminus-(Aself-eye(3*S.npts)),'fro')/norm(Aself,'fro')<5e-14)
assert(relative_error(self_minus_fmm,self_fmm-density)<5e-14)
assert(relative_error(self_minus_fmm,self_minus)<5e-8)

%% Source-to-target: rectangular matrix, no identity, FMM and dense agree.
[cross_dense,Across] = nrccie_apply( ...
    S,density,T,Ccross,dense_options);
assert(isequal(size(cross_dense),[3,T.npts]))
assert(isequal(size(Across),[3*T.npts,3*S.npts]))
assert(relative_error(cross_dense(:),Across*density(:))<5e-14)
cross_vector = nrccie_apply(S,density(:),T,Ccross,negative_dense);
assert(isequal(size(cross_vector),[3*T.npts,1]))
assert(relative_error(cross_vector,cross_dense(:))<5e-14)

[cross_fmm,Across_fmm] = nrccie_apply( ...
    S,density,T,Ccross,negative_fmm);
assert(isempty(Across_fmm))
assert(relative_error(cross_fmm,cross_dense)<5e-8)

%% Zero source weights isolate the self jump from the principal value.
zero_source.r = S.r;
zero_source.n = S.n;
zero_source.npts = S.npts;
zero_source.dru = S.dru;
zero_source.drv = S.drv;
zero_source.wts = zeros(1,S.npts);
[jump_only,Ajump] = nrccie_apply( ...
    zero_source,density,[],dense_options);
assert(relative_error(jump_only,0.5*density)<5e-14)
assert(norm(Ajump-0.5*eye(3*S.npts),'fro')<5e-14)
[cross_zero,Across_zero] = nrccie_apply( ...
    zero_source,density,T,[],dense_options);
assert(norm(cross_zero,'fro')==0)
assert(norm(Across_zero,'fro')==0)

fprintf(['NRCCIE_APPLY_PASS: self FMM/dense %.3e, ' ...
    'cross FMM/dense %.3e\n'], ...
    relative_error(self_fmm,self_dense), ...
    relative_error(cross_fmm,cross_dense))


function [potential,gx,gy,gz] = apply_full_correction(C,density)
potential = (C{1}*density.').';
gx = (C{2}*density.').';
gy = (C{3}*density.').';
gz = (C{4}*density.').';
end


function error = relative_error(a,b)
scale = max(norm(b(:)),realmin);
error = norm(a(:)-b(:))/scale;
end
