function rotations = nrccie_deterministic_rotations(count,phase,fourfold)
% Deterministic incident directions with two transverse polarizations.
% No RNG or random projection. A phase offset separates training/check sets.
% With fourfold=true, count/4 seed frames generate complete C4 orbits;
% count is still the TOTAL number of propagation directions.
if nargin<2, phase = 0; end
if nargin<3, fourfold = false; end
validateattributes(count,{'numeric'},{'scalar','integer','positive','finite'})
validateattributes(phase,{'numeric'},{'scalar','real','finite','>=',0,'<',1})
validateattributes(fourfold,{'logical'},{'scalar'})
if fourfold
    assert(mod(count,4)==0,'C4 training needs a multiple of four directions.')
    count = count/4;
end
j = (0:count-1)+0.5;
z = 1-2*j/count;
theta = pi*(3-sqrt(5))*j+2*pi*phase;
radial = sqrt(max(0,1-z.^2));
d = [radial.*cos(theta);radial.*sin(theta);z];
e1 = [-sin(theta);cos(theta);zeros(1,count)];
e2 = cross(d,e1,1);
rotations = zeros(3,3,count);
rotations(:,1,:) = reshape(e1,3,1,count);
rotations(:,2,:) = reshape(e2,3,1,count);
rotations(:,3,:) = reshape(d,3,1,count);
if fourfold
    Q = [0,1,0;-1,0,0;0,0,1];
    base = rotations;
    rotations = zeros(3,3,4*count);
    for p = 0:3
        rotations(:,:,p*count+(1:count)) = pagemtimes(Q^p,base);
    end
end
end
