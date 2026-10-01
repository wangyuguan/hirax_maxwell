function rotations = sample_so3(nsample)
%SAMPLE_SO3 Haar-uniform rotations, returned as a 3-by-3-by-nsample array.
% Uses the caller's RNG state. Gaussian QR needs positive diagonal signs.
rotations = zeros(3,3,nsample);
for j = 1:nsample
    [Q,R] = qr(randn(3));
    Q = Q*diag(sign(diag(R)));
    Q(:,3) = sign(det(Q))*Q(:,3);
    rotations(:,:,j) = Q;
end
end
