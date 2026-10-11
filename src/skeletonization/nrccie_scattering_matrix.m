function [Scat,info] = nrccie_scattering_matrix(F,Afun,tol,max_refinement, ...
        blocks_in,Pmerge_in,blocks_out,Pmerge_out,batch_size)
%NRCCIE_SCATTERING_MATRIX Form Pout*A^{-1}*Pin with direct RS-S solves.
%
% Input arguments:
%   F - native direct RS-S factor of the complete corrected A.
%   Afun - application of the original corrected A for residual refinement.
%   tol, max_refinement - requested global residual and correction limit.
%   blocks_in, Pmerge_in - multilevel incoming interpolation maps.
%   blocks_out, Pmerge_out - multilevel outgoing equivalent-density maps.
%   batch_size - number of incoming columns solved together.
%
% Output arguments:
%   Scat - 3*Nout-by-3*Nin matrix Pout*A^{-1}*Pin.
%   info - per-batch direct-solver diagnostics and assembly time.
%
% With Jin/Jout from the skeleton functions, qout = Scat*b(Jin,:) and
% fields = K(:,Jout)*qout. K includes source quadrature weights; qout is
% an equivalent density, not the restriction of the solved density.
% Full interpolation matrices and the system inverse are never formed.
validateattributes(batch_size,{'numeric'}, ...
    {'scalar','finite','integer','positive'})
nin = size(Pmerge_in,2);
nout = size(Pmerge_out,1);
Scat = complex(zeros(nout,nin));
info = struct('batches',{{}},'columns',{{}},'seconds',0);
timer = tic;
batch = 0;
for first = 1:batch_size:nin
    columns = first:min(first+batch_size-1,nin);
    rhs = complex(zeros(F.N,numel(columns)));
    for part = 1:numel(blocks_in)
        block = blocks_in(part);
        dofs = reshape(3*block.nodes+(-2:0).',1,[]);
        rhs(dofs,:) = block.P*Pmerge_in(block.candidate_dofs,columns);
    end
    batch = batch+1;
    [density,batch_info] = nrccie_direct_solve( ...
        F,rhs,Afun,tol,max_refinement);
    info.batches{batch} = batch_info;
    info.columns{batch} = columns;
    for part = 1:numel(blocks_out)
        block = blocks_out(part);
        dofs = reshape(3*block.nodes+(-2:0).',1,[]);
        Scat(:,columns) = Scat(:,columns)+ ...
            Pmerge_out(:,block.candidate_dofs)*(block.P*density(dofs,:));
    end
    fprintf(['Scat columns %d:%d/%d: direct RS-S, refinement %d/%d, ' ...
        'max A residual %.3e, elapsed %.2f s.\n'], ...
        columns(1),columns(end),nin,min(batch_info.refinement_steps), ...
        max(batch_info.refinement_steps),max(batch_info.true_relres),toc(timer))
end
info.seconds = toc(timer);
end
