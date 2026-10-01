% Singular values of the matrix saved by run_two_leaf_system_matrix.m.
% Default: full two-leaf system. Assemble it first with that script.
% NRCCIE_TWO_LEAF_ROWS_ONLY=1 instead analyzes the saved sampled rows;
% their spectrum/rank/condition number do not describe the full system.
% NRCCIE_TWO_LEAF_SVD_VECTORS=1 also computes and saves U,Sigma,V.
% A full run is expensive: A alone occupies 4.67 GiB at order 1, and
% computing vectors adds two similarly sized matrices plus workspace.
% These are the raw coefficient-matrix singular values in the Euclidean
% DOF norm; no quadrature-weighted L2 rescaling is applied.
clear
root = fileparts(fileparts(mfilename('fullpath')));
leaf_ids = [1,2];
norder = 1;
rows_only = strcmp(getenv('NRCCIE_TWO_LEAF_ROWS_ONLY'),'1');
compute_vectors = strcmp(getenv('NRCCIE_TWO_LEAF_SVD_VECTORS'),'1');
suffix = '';
if rows_only, suffix = '_rows'; end
stem = sprintf('two_leaf_system_%d_%d_order%d%s', ...
    leaf_ids(1),leaf_ids(2),norder,suffix);
input_file = fullfile(root,'data',[stem,'.mat']);
assert(isfile(input_file), ...
    'Matrix file missing. Run run_two_leaf_system_matrix.m first.')
% Do not load surfer objects or start the quadrature/FMM dependencies.
load(input_file,'assembly','validation','matrix_rows','settings')
assert(assembly.complete && validation.passed, ...
    'Matrix assembly/verification has not finished successfully.')
assert(assembly.rows_only==rows_only,'Unexpected matrix provenance.')
if ~rows_only
    assert(isequal(assembly.matrix_size,assembly.full_matrix_size) ...
        && isequal(matrix_rows,1:assembly.full_matrix_size(1)), ...
        'The saved matrix is not the complete system.')
end
load(input_file,'A')
assert(isequal(size(A),assembly.matrix_size) && all(isfinite(A),'all'))
matrix_size = size(A);
scope = 'full system';
if rows_only, scope = 'sampled rows only'; end
fprintf('SVD of %s: %d x %d complex matrix.\n',scope,matrix_size)
timer = tic;
if compute_vectors
    [U,Sigma,V] = svd(A,'econ');
    singular_values = diag(Sigma);
else
    singular_values = svd(A);
end
svd_seconds = toc(timer);
relative_singular_values = singular_values/singular_values(1);
fprintf('SVD time %.2f s; sigma_max %.8e; sigma_min %.8e.\n', ...
    svd_seconds,singular_values(1),singular_values(end))
if rows_only
    fprintf('These are singular values of selected rows, not the full system.\n')
end
output_file = fullfile(root,'data',[stem,'_svd.mat']);
save(output_file,'singular_values','relative_singular_values', ...
    'svd_seconds','matrix_size','matrix_rows','settings','assembly', ...
    'input_file','scope','compute_vectors','-v7.3')
if compute_vectors
    save(output_file,'U','Sigma','V','-append')
end
fig = figure('Visible','off','Color','w');
semilogy(1:numel(singular_values),singular_values,'o-','MarkerSize',3)
xlabel('Index'); ylabel('Singular value'); grid on
title(sprintf('Two-leaf NRCCIE: %s (%d x %d)',scope,matrix_size))
plot_file = fullfile(root,'data',[stem,'_svd.png']);
exportgraphics(fig,plot_file,'Resolution',180)
close(fig)
fprintf('Saved %s\nSaved %s\n',output_file,plot_file)
