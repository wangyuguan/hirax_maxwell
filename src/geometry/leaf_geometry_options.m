function o = leaf_geometry_options(o)
% Shared mesh defaults for four copies of the original closed leaf.
if nargin<1 || isempty(o), o = struct(); end
if ~isfield(o,'norder'), o.norder = 6; end
if ~isfield(o,'rim_width'), o.rim_width = 0.028128271246*69.25; end
if ~isfield(o,'cap_collar_width'), o.cap_collar_width = o.rim_width; end
if ~isfield(o,'cap_mesh_spacing'), o.cap_mesh_spacing = 0.2*69.25; end
if ~isfield(o,'cap_mesh_side_start'), o.cap_mesh_side_start = 0.62; end
if ~isfield(o,'outline_refinement'), o.outline_refinement = 1; end
if ~isfield(o,'wall_profile_refinement'), o.wall_profile_refinement = 3; end
if ~isfield(o,'surface_gap'), o.surface_gap = 4; end
end
