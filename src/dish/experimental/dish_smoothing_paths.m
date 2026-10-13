function [mesh_dir,smoother_root] = dish_smoothing_paths(stem)
%DISH_SMOOTHING_PATHS Locate experimental inputs independently of the cwd.

code_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(fileparts(code_dir)));
dependency_root = getenv('NRCCIE_DEPENDENCY_ROOT');
if isempty(dependency_root)
    dependency_root = fileparts(project_root);
end
smoother_root = getenv('SURFSMOOTH3D_ROOT');
if isempty(smoother_root)
    smoother_root = fullfile(dependency_root,'SurfSmooth3D');
end
mesh_root = getenv('DISH_MESH_ROOT');
if isempty(mesh_root)
    mesh_root = fullfile(project_root,'data','dish');
end
mesh_dir = fullfile(mesh_root,stem);
end
