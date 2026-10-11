function [X,info] = nrccie_direct_solve(F,B,Afun,tol,max_refinement)
%NRCCIE_DIRECT_SOLVE Solve A*X=B with a direct RS-S factor of A.
% The native backsolve is followed, when necessary, by residual corrections
% X <- X + F^{-1}(B-A*X). This is direct refinement, not GMRES.

if nargin<5, max_refinement = 3; end
validateattributes(B,{'numeric'},{'2d','nrows',F.N,'finite'})
validateattributes(tol,{'numeric'}, ...
    {'scalar','real','finite','positive','<',1})
validateattributes(max_refinement,{'numeric'}, ...
    {'scalar','real','finite','integer','nonnegative'})
assert(isa(Afun,'function_handle'),'Afun must apply the original operator A.')

total_timer = tic;
timer = tic;
X = srskelf_sv_nn(F,B);
backsolve_seconds = toc(timer);
[R,relres,residual_seconds] = residual(Afun,X,B);
initial_relres = relres;
steps = zeros(1,size(B,2));
refinement_seconds = 0;

for step = 1:max_refinement
    active = find(~isfinite(relres) | relres>tol);
    if isempty(active), break; end
    timer = tic;
    X(:,active) = X(:,active)+srskelf_sv_nn(F,R(:,active));
    refinement_seconds = refinement_seconds+toc(timer);
    steps(active) = step;
    [R,relres,part_seconds] = residual(Afun,X,B);
    residual_seconds = residual_seconds+part_seconds;
end

failed = find(~isfinite(relres) | relres>tol);
if ~isempty(failed)
    error('nrccie_direct_solve:NoConvergence', ...
        ['Direct RS-S refinement failed for RHS %d after %d corrections: ' ...
         'relative residual %.3e > %.3e.'], ...
        failed(1),steps(failed(1)),relres(failed(1)),tol)
end
info = struct('method','direct-a-rss','initial_true_relres',initial_relres, ...
    'true_relres',relres,'refinement_steps',steps, ...
    'backsolve_seconds',backsolve_seconds, ...
    'refinement_seconds',refinement_seconds, ...
    'residual_seconds',residual_seconds,'seconds',toc(total_timer));
end


function [R,relres,seconds] = residual(Afun,X,B)
timer = tic;
R = B-Afun(X);
norm_b = vecnorm(B,2,1);
relres = vecnorm(R,2,1)./max(norm_b,realmin);
zero = norm_b==0;
relres(zero) = vecnorm(R(:,zero),2,1);
seconds = toc(timer);
end
