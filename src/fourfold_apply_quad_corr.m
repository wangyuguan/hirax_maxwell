function [potential,gradient_x,gradient_y,gradient_z] = ...
    fourfold_apply_quad_corr(corrections,density)
%FOURFOLD_APPLY_QUAD_CORR Apply the block-circulant quadrature correction.
%
% DENSITY is (ndensity, 4*component_npts) in the merge order used to build
% CORRECTIONS, and the four outputs have the same size. The result equals
% multiplying by the assembled matrices of FOURFOLD_EXPAND_QUAD_CORR, but
% only the single stored block row is touched.
%
% The block row is held with S1 as its target. For any other target the
% same blocks apply once the accumulated gradient is rotated out of the S1
% frame, which is what the trailing rotation does.

component_npts = corrections.component_npts;
cycle_to_merge = corrections.cycle_to_merge;
Q = corrections.rotation;

ndensity = size(density,1);

potential = complex(zeros(ndensity,4*component_npts));
gradient_x = complex(zeros(ndensity,4*component_npts));
gradient_y = complex(zeros(ndensity,4*component_npts));
gradient_z = complex(zeros(ndensity,4*component_npts));

for target_cycle = 1:4
    accumulated_potential = complex(zeros(component_npts,ndensity));
    accumulated_x = complex(zeros(component_npts,ndensity));
    accumulated_y = complex(zeros(component_npts,ndensity));
    accumulated_z = complex(zeros(component_npts,ndensity));

    for source_cycle = 1:4
        relative_cycle = mod(source_cycle-target_cycle,4)+1;
        if ~corrections.nonempty(relative_cycle)
            continue
        end
        source_columns = component_columns( ...
            cycle_to_merge(source_cycle),component_npts);
        source_density = density(:,source_columns).';

        if ~isempty(corrections.reflection)
            B = {corrections.slp{relative_cycle},corrections.gx{relative_cycle}, ...
                corrections.gy{relative_cycle},corrections.gz{relative_cycle}};
            [p,x,y,z] = apply_midplane_quad_corr(B,source_density, ...
                corrections.reflection,corrections.reflection);
            accumulated_potential = accumulated_potential+p;
            accumulated_x = accumulated_x+x;
            accumulated_y = accumulated_y+y;
            accumulated_z = accumulated_z+z;
            continue
        end

        accumulated_potential = accumulated_potential+ ...
            corrections.slp{relative_cycle}*source_density;
        accumulated_x = accumulated_x+ ...
            corrections.gx{relative_cycle}*source_density;
        accumulated_y = accumulated_y+ ...
            corrections.gy{relative_cycle}*source_density;
        accumulated_z = accumulated_z+ ...
            corrections.gz{relative_cycle}*source_density;
    end

    % The stored gradient blocks are expressed in the S1 frame; rotate
    % them into the frame of this target surface.
    Qtarget = Q^(target_cycle-1);
    target_columns = component_columns( ...
        cycle_to_merge(target_cycle),component_npts);

    potential(:,target_columns) = accumulated_potential.';
    gradient_x(:,target_columns) = (Qtarget(1,1)*accumulated_x+ ...
        Qtarget(1,2)*accumulated_y+Qtarget(1,3)*accumulated_z).';
    gradient_y(:,target_columns) = (Qtarget(2,1)*accumulated_x+ ...
        Qtarget(2,2)*accumulated_y+Qtarget(2,3)*accumulated_z).';
    gradient_z(:,target_columns) = (Qtarget(3,1)*accumulated_x+ ...
        Qtarget(3,2)*accumulated_y+Qtarget(3,3)*accumulated_z).';
end
end


function columns = component_columns(merge_index,component_npts)
columns = (merge_index-1)*component_npts+(1:component_npts);
end
