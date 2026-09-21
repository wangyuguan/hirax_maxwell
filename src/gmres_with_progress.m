function [x,flag,relres,iter,resvec] = gmres_with_progress(afun,b,restart,tol,maxit)
% GMRES_WITH_PROGRESS Run MATLAB GMRES with timed operator-call logging.
% AFUN is a function handle. Arguments and outputs match the five-argument
% GMRES call. Operator evaluations include residual checks, so their count
% is not an iteration count. No extra operator evaluations are performed.

started = tic;
calls = 0;
fprintf('GMRES started: tolerance %.3e, restart %s, maxit %d\n', ...
    tol,mat2str(restart),maxit);
[x,flag,relres,iter,resvec] = gmres(@logged_matvec,b,restart,tol,maxit);
fprintf(['GMRES finished: flag %d, iteration %s, relative residual %.3e, ' ...
    '%d matvecs, elapsed %.1f s\n'], ...
    flag,mat2str(iter),relres,calls,toc(started));

    function y = logged_matvec(v)
        calls = calls+1;
        fprintf('GMRES matvec %d started; elapsed %.1f s\n',calls,toc(started));
        call_started = tic;
        y = afun(v);
        fprintf('GMRES matvec %d completed in %.1f s; elapsed %.1f s\n', ...
            calls,toc(call_started),toc(started));
    end
end
