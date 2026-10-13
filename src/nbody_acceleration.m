function a = nbody_acceleration(r, GM)
%NBODY_ACCELERATION  Mutual Newtonian accelerations of N massive bodies.
%
%   a = NBODY_ACCELERATION(r, GM) with positions r (3xN) [m] and
%   gravitational parameters GM (1xN) [m^3 s^-2] returns the accelerations
%   a (3xN) [m s^-2]:
%
%       a_i = - sum_{j ~= i} GM_j (r_i - r_j) / |r_i - r_j|^3
%
%   Direct O(N^2) summation, vectorised over all pairs.

n  = size(r, 2);
dx = r(1, :).' - r(1, :);          % dx(i,j) = x_i - x_j
dy = r(2, :).' - r(2, :);
dz = r(3, :).' - r(3, :);
inv_d3 = (dx.^2 + dy.^2 + dz.^2) .^ (-1.5);
inv_d3(1:n+1:end) = 0;              % no self-interaction
w = GM(:).' .* inv_d3;              % w(i,j) = GM_j / d_ij^3
a = -[sum(w .* dx, 2), sum(w .* dy, 2), sum(w .* dz, 2)].';
end
