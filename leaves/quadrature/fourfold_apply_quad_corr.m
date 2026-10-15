function [potential,gradient_x,gradient_y,gradient_z] = ...
    fourfold_apply_quad_corr(corrections,density)
% FOURFOLD_APPLY_QUAD_CORR Apply the C4 quadrature correction.
%
%  Syntax
%   [p,gx,gy,gz] = fourfold_apply_quad_corr(C,density)
%
%  Input arguments:
%    * C: correction struct returned by fourfold_quad_corr_mats
%    * density: ndensity-by-C.npts array in the merged surface order
%
%  Output arguments:
%    * p,gx,gy,gz: potential and Cartesian gradient corrections, each with
%        the same size as density
%
%  Reuse the stored target row and rotate the gradient to each target frame.

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
        source_columns = (cycle_to_merge(source_cycle)-1)*component_npts+ ...
            (1:component_npts);
        source_density = density(:,source_columns).';

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
      target_columns = (cycle_to_merge(target_cycle)-1)*component_npts+ ...
          (1:component_npts);

      potential(:,target_columns) = accumulated_potential.';
      gradient_x(:,target_columns) = (Qtarget(1,1)*accumulated_x+ ...
          Qtarget(1,2)*accumulated_y+Qtarget(1,3)*accumulated_z).';
      gradient_y(:,target_columns) = (Qtarget(2,1)*accumulated_x+ ...
          Qtarget(2,2)*accumulated_y+Qtarget(2,3)*accumulated_z).';
      gradient_z(:,target_columns) = (Qtarget(3,1)*accumulated_x+ ...
          Qtarget(3,2)*accumulated_y+Qtarget(3,3)*accumulated_z).';
    end
end
