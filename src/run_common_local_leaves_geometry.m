function [S,parts] = run_common_local_leaves_geometry(norder,opts)
% Save one common mesh, four local meshes, and their rigid transforms.

if nargin<1 || isempty(norder), norder=6; end
if nargin<2, opts=struct(); end
root=fileparts(fileparts(mfilename('fullpath')));
load_dependencies(root);
addpath(fullfile(root,'src'));
opts.norder=norder;
thickness=.6925;
if isfield(opts,'thickness'), thickness=opts.thickness; end
output_file=fullfile(root,'data',sprintf('common_local_geometry_order%02d.mat',norder));
if isfield(opts,'output_file'), output_file=opts.output_file; end
[S,parts]=hirax_common_local_leaves_surfer(thickness,opts);
geometry=struct('common',parts.local_common_surfer, ...
    'local_tips',{parts.local_unique_surfers}, ...
    'rotation',parts.rotation,'shift',parts.shift);
metadata=rmfield(parts,{'local_common_surfer','local_unique_surfers','leaves'});
output_dir=fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir), mkdir(output_dir); end
save(output_file,'geometry','metadata','-v7.3');
end

function load_dependencies(root)
% Isolate variables introduced by the dependency startup scripts.
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'));
run(fullfile(root,'..','chunkie','startup.m'));
end
