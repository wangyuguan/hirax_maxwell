function o = leaf_geometry_options(o)
% Shared defaults for the straight-interface four-tip geometry.
if nargin<1 || isempty(o), o = struct(); end
if isfield(o,'extend_common_core') && o.extend_common_core
    error('Curved common-core extension is retired; use the straight interface.');
end
if isfield(o,'common_core_clearance') || isfield(o,'common_interface_tolerance')
    error('Remove retired curved-interface options; the current partition is straight.');
end
if ~isfield(o,'extend_common_outline'), o.extend_common_outline = true; end
if ~isfield(o,'norder'), o.norder = 6; end
if ~isfield(o,'common_rim_width'), o.common_rim_width = .15; end
if ~isfield(o,'common_collar_width'), o.common_collar_width = .6; end
if ~isfield(o,'tip_rim_width'), o.tip_rim_width = []; end
if ~isfield(o,'tip_collar_width'), o.tip_collar_width = []; end
if ~isfield(o,'width_transition_y'), o.width_transition_y = [-55,-40]; end
if ~isfield(o,'cap_mesh_spacing'), o.cap_mesh_spacing = 6; end
if ~isfield(o,'local_cap_mesh_spacing'), o.local_cap_mesh_spacing = 2; end
if ~isfield(o,'interface_spacing'), o.interface_spacing = 4; end
if ~isfield(o,'wall_profile_refinement'), o.wall_profile_refinement = 3; end
if ~isfield(o,'surface_gap'), o.surface_gap = 4; end
if ~isfield(o,'outline_geometry_tolerance'), o.outline_geometry_tolerance = 1e-6; end
validateattributes(o.interface_spacing,{'numeric'},{'scalar','real','finite','positive'});
end
