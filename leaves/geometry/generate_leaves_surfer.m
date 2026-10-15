function [S,parts] = generate_leaves_surfer(thickness,opts)
% GENERATE_LEAVES_SURFER Construct four identical closed leaves.
%
%  Syntax
%   [S,parts] = generate_leaves_surfer(thickness,opts)
%
%  Input arguments:
%    * thickness: leaf thickness in mm
%    * opts: mesh options for leaf_plate_surfer, with surface_gap in mm
%
%  Output arguments:
%    * S: merged surfer in top-left, top-right, bottom-left, bottom-right order
%    * parts: individual surfaces, node/patch indices, rotations, and shifts
%
%  The master leaf is meshed once and copied by rigid motions.

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
