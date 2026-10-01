function Scat = nrccie_scattering_matrix(F,blocks_in,Pmerge_in,blocks_out,Pmerge_out,batch_size)
%NRCCIE_SCATTERING_MATRIX Map incoming skeleton RHS to outgoing density.
% F is returned by nrccie_rskelf; blocks and Pmerge are returned by the
% corresponding hierarchical_nrccie_incoming/outgoing_skeleton functions.
% Scat = Pout * A^(-1) * Pin, with size 3*Nout-by-3*Nin.
% For the Jin/Jout ordering from those skeleton functions,
%   qout = Scat*b(Jin,:);
%   fields = K(:,Jout)*qout;
% K maps the full source density to fields, including source weights.
% qout is an equivalent density, not a restriction of the solved density.
% Only batch_size columns of Pin and their solutions are expanded at once;
% neither the full interpolation matrices nor the inverse of A are formed.
nin = size(Pmerge_in,2);
nout = size(Pmerge_out,1);
Scat = complex(zeros(nout,nin));
for first = 1:batch_size:nin
    columns = first:min(first+batch_size-1,nin);
    rhs = complex(zeros(F.N,numel(columns)));
    for part = 1:numel(blocks_in)
        block = blocks_in(part);
        dofs = reshape(3*block.nodes+(-2:0).',1,[]);
        rhs(dofs,:) = block.P*Pmerge_in(block.candidate_dofs,columns);
    end
    density = rskelf_sv(F,rhs);
    for part = 1:numel(blocks_out)
        block = blocks_out(part);
        dofs = reshape(3*block.nodes+(-2:0).',1,[]);
        Scat(:,columns) = Scat(:,columns)+ ...
            Pmerge_out(:,block.candidate_dofs)*(block.P*density(dofs,:));
    end
end
end
