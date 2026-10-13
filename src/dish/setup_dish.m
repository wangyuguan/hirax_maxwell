function setup_dish()
%SETUP_DISH Add dish geometry tools and fmm3dbie to the MATLAB path.
%   Optional meshing and smoothing dependencies are loaded by their drivers.

dish_root = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(dish_root));
dependency_root = getenv('NRCCIE_DEPENDENCY_ROOT');
if isempty(dependency_root)
    dependency_root = fileparts(project_root);
end

addpath(dish_root);
addpath(fullfile(dish_root,'mesh'));
if isempty(which('surfer')) || isempty(which('koorn.rv_nodes'))
    fmm3dbie_path = fullfile(dependency_root,'fmm3dbie-hirax-dev','matlab');
    assert(isfolder(fmm3dbie_path), ...
        'Missing fmm3dbie MATLAB folder %s; set NRCCIE_DEPENDENCY_ROOT.', ...
        fmm3dbie_path);
    addpath(fmm3dbie_path);
    addpath(fullfile(fmm3dbie_path,'src'));
end
end
