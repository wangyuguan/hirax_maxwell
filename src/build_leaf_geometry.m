function [geometry,metadata,S,parts] = build_leaf_geometry(norder,opts)
% Build current geometry rather than silently loading a historical MAT mesh.
if nargin<1 || isempty(norder), norder = 6; end
if nargin<2, opts = struct(); end
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'..','fmm3dbie-hirax-dev','matlab','startup.m'))
run(fullfile(root,'..','chunkie','startup.m'))
opts.norder = norder;
thickness = .6925;
if isfield(opts,'thickness'), thickness = opts.thickness; end
[S,parts] = generate_leaves_surfer(thickness,opts);
geometry.common = parts.local_common_surfer;
geometry.local_tips = parts.local_unique_surfers;
geometry.rotation = parts.rotation;
geometry.shift = parts.shift;
metadata = rmfield(parts,{'local_common_surfer','local_unique_surfers','leaves'});
end
