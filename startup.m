% Dependencies and HIRAX MATLAB routines.
hirax_root = fileparts(mfilename('fullpath'));
dependency_root = getenv('NRCCIE_DEPENDENCY_ROOT');
if isempty(dependency_root), dependency_root = fileparts(hirax_root); end
run(fullfile(dependency_root,'fmm3dbie-hirax-dev','matlab','startup.m'))
% Only Chunkie's geometry routines are needed; do not initialize/build FMM2D.
addpath(fullfile(dependency_root,'chunkie','chunkie'))
addpath(fullfile(dependency_root,'fmm3dbie-hirax-dev','FMM3D','matlab'))
addpath(genpath(fullfile(dependency_root,'FLAM')))
addpath(genpath(fullfile(hirax_root,'src')))
addpath(fullfile(dependency_root,'strong-skel'))
addpath(fullfile(hirax_root,'test'))
addpath(fullfile(hirax_root,'test','ellipsoid'))
clear hirax_root dependency_root pth dir
