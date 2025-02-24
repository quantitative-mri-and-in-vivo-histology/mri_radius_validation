function angles = compute_angles(a,b)
% Computes pairwise angles between 3D vectors in matrices a and b.
%
% INPUTS:
%   a - (Nx3 matrix) Set of N 3D vectors.
%   b - (Mx3 matrix) Set of M 3D vectors.
%
% OUTPUT:
%   angles - (NxM matrix) Matrix of angles (radians) between each pair of vectors.

    [A,B] = meshgrid(1:size(a,1), 1:size(b,1));
    A_lin = A(:);
    B_lin = B(:);
    N1 = a(A_lin,:);
    N2 = b(B_lin,:);
    N1dotN2 = N1(:, 1) .* N2(:, 1) + N1(:, 2) .* N2(:, 2) + N1(:, 3) .* N2(:, 3);
    N1xN2   = [(N1(:, 2) .* N2(:, 3) - N1(:, 3) .* N2(:, 2)), ...
         (N1(:, 3) .* N2(:, 1) - N1(:, 1) .* N2(:, 3)), ...
         (N1(:, 1) .* N2(:, 2) - N1(:, 2) .* N2(:, 1))];
    angles   = atan2(sqrt(sum(N1xN2 .* N1xN2, 2)), N1dotN2);
    angles = reshape(angles, size(b,1), size(a,1))';
    angles = reshape(angles, size(a,1), size(b,1));
end
