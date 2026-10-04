function [S,parts] = generate_leaves_surfer(thickness,opts)
% Four identical closed leaves in TL/TR/BL/BR order; lengths in mm.
% Mesh the original leaf once and rotate/translate copies of that mesh.
% S is merged in leaf order; parts.node_ids/patch_ids map into S.
if nargin<2 || isempty(opts), opts = struct(); end
opts = leaf_geometry_options(opts);
validateattributes(opts.norder,{'numeric'},{'scalar','integer','>=',1})
[base,base_parts] = leaf_plate_surfer(thickness,opts);

parts.base_surfer = base;
parts.base_parts = base_parts;
parts.geometry_revision = leaf_geometry_revision();
parts.thickness = thickness;
parts.options = base_parts.options;
parts.leaf_ids = 1:4;
parts.leaves = cell(1,4);
parts.node_ids = cell(1,4);
parts.patch_ids = cell(1,4);
parts.rotation = zeros(3,3,4);
parts.shift = zeros(3,4);

% The original long flanks have a 2.83*sqrt(2) mm separation.
center = base_parts.design.leaf_center_raw;
center(2) = center(2)+(opts.surface_gap-2.83*sqrt(2))/sqrt(2);
angles = [pi/4,-pi/4,3*pi/4,-3*pi/4];
for leaf = 1:4
    a = angles(leaf);
    R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    shift = R*[center;0];
    parts.leaves{leaf} = affine_transf(base,R,shift);
    parts.rotation(:,:,leaf) = R;
    parts.shift(:,leaf) = shift;
    parts.node_ids{leaf} = (leaf-1)*base.npts+(1:base.npts);
    parts.patch_ids{leaf} = (leaf-1)*base.npatches+(1:base.npatches);
end
S = merge([parts.leaves{:}]);
end
