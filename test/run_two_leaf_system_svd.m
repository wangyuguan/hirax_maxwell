% SERVER RUN: assemble the complete two-leaf system, then compute its SVD.
% From the repository root:
%   matlab -batch "run('test/run_two_leaf_system_svd.m')"
%
% The assembly needs sibling directories ../fmm3dbie-hirax-dev and
% ../chunkie, with their MATLAB paths and server-native MEX dependencies.
% Geometry/physics settings are in run_two_leaf_system_matrix.m; defaults
% use original leaves 1 and 2, order 1, and 17712 scalar unknowns.
% A is streamed to data/two_leaf_system_1_2_order1.mat (about 4.67 GiB).
% SVD loads A into RAM and needs additional workspace. By default it saves
% all singular values and a spectrum plot to data/*_svd.mat and *_svd.png.
% To also save U,Sigma,V, set NRCCIE_TWO_LEAF_SVD_VECTORS=1 before running.
% To repeat only the SVD, run test/run_two_leaf_svd.m separately.
%
% This entry point always runs the FULL system, even after a local smoke test.
setenv('NRCCIE_TWO_LEAF_ROWS_ONLY','0')
run(fullfile(fileparts(mfilename('fullpath')),'run_two_leaf_system_matrix.m'))
run(fullfile(fileparts(mfilename('fullpath')),'run_two_leaf_svd.m'))
