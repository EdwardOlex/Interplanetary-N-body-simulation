function a = test_particle_acceleration(rp, rb, GM)
%TEST_PARTICLE_ACCELERATION  Acceleration of massless particles in the field of massive bodies.
%
%   a = TEST_PARTICLE_ACCELERATION(rp, rb, GM) returns the accelerations
%   a (3xM) [m s^-2] of M test particles at positions rp (3xM) produced by
%   N massive bodies at positions rb (3xN) with gravitational parameters
%   GM (1xN). The particles do not act on the massive bodies (restricted
%   N-body problem), which keeps the cost linear in the number of particles.

dx = rp(1, :).' - rb(1, :);         % (M x N)
dy = rp(2, :).' - rb(2, :);
dz = rp(3, :).' - rb(3, :);
w  = GM(:).' .* (dx.^2 + dy.^2 + dz.^2) .^ (-1.5);
a  = -[sum(w .* dx, 2), sum(w .* dy, 2), sum(w .* dz, 2)].';
end
