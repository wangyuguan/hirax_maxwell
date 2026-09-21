function [S,parts] = run_common_local_leaves_geometry(norder,opts)
% Save one common mesh, four local meshes, and their rigid transforms.

if nargin<1 || isempty(norder), norder = 6; end
if nargin<2, opts = struct(); end
root = fileparts(fileparts(mfilename('fullpath')));
output_file = fullfile(root,'data',sprintf('common_local_geometry_order%02d.mat',norder));
if isfield(opts,'output_file'), output_file = opts.output_file; end
[geometry,metadata,S,parts] = build_leaf_geometry(norder,opts);
output_dir = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir), mkdir(output_dir); end
save(output_file,'geometry','metadata','-v7.3');
end
