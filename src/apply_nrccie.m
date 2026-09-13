function y = apply_nrccie(x,operator)
% Unknowns [ju;jv;rho] are interleaved by node. The four scalar kernels
% interpolate Cartesian current components and charge on each patch.

npts = operator.npts;
density_components = reshape(x,3,npts);
surface_current = operator.ru.*density_components(1,:)+ ...
    operator.rv.*density_components(2,:);
surface_charge = density_components(3,:);
density = [surface_current;surface_charge];

source = struct();
source.sources = operator.r;
source.nd = 4;
source.charges = density.*operator.wts;

fmm_output = hfmm3d(operator.eps_fmm,operator.zk,source,2);
potential = reshape(fmm_output.pot,4,npts);
gradient = reshape(fmm_output.grad,4,3,npts);

[potential_correction,gradient_x,gradient_y,gradient_z] = ...
    operator.apply_corrections(density);
potential = potential+potential_correction;
gradient(:,1,:) = gradient(:,1,:)+reshape(gradient_x,4,1,npts);
gradient(:,2,:) = gradient(:,2,:)+reshape(gradient_y,4,1,npts);
gradient(:,3,:) = gradient(:,3,:)+reshape(gradient_z,4,1,npts);

slp_current = potential(1:3,:);
slp_charge = potential(4,:);
gradient_slp_charge = reshape(gradient(4,:,:),3,npts);

curl_slp_current = complex(zeros(3,npts));
curl_slp_current(1,:) = reshape( ...
    gradient(3,2,:)-gradient(2,3,:),1,npts);
curl_slp_current(2,:) = reshape( ...
    gradient(1,3,:)-gradient(3,1,:),1,npts);
curl_slp_current(3,:) = reshape( ...
    gradient(2,1,:)-gradient(1,2,:),1,npts);
divergence_slp_current = reshape( ...
    gradient(1,1,:)+gradient(2,2,:)+gradient(3,3,:),1,npts);

electric_field = 1i*operator.zk*slp_current-gradient_slp_charge;
nxh = cross(operator.n,curl_slp_current,1);
normal_electric_field = sum(operator.n.*electric_field,1);
nxnxe = operator.n.*normal_electric_field-electric_field;

principal_value = complex(zeros(3,npts));
tangent_equation = -nxh+operator.alpha*nxnxe;
principal_value(1,:) = sum(operator.ru.*tangent_equation,1);
principal_value(2,:) = sum(operator.rv.*tangent_equation,1);
principal_value(3,:) = -normal_electric_field+operator.alpha*( ...
    divergence_slp_current-1i*operator.zk*slp_charge);

y = 0.5*density_components+principal_value;
y = y(:);
end
