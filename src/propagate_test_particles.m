function part = propagate_test_particles(traj, rp0, vp0)
%PROPAGATE_TEST_PARTICLES  Propagate massless bodies (e.g. asteroids) through a precomputed ephemeris.
%
%   part = PROPAGATE_TEST_PARTICLES(traj, rp0, vp0) integrates M test
%   particles with initial positions RP0 (3xM) [m] and velocities VP0 (3xM)
%   [m/s] in the gravitational field of the massive bodies stored in TRAJ
%   (output of PROPAGATE_NBODY), on the same time grid and with the same
%   explicit midpoint scheme.
%
%   The massive bodies are propagated first and their ephemeris is reused
%   for every particle: particles feel the planets but do not perturb them.
%   Planet positions at the half step are obtained by linear interpolation
%   between consecutive ephemeris epochs.
%
%   Output struct PART
%     t   1xT epochs [s]
%     r   3xMxT positions [m]
%     v   3xM final velocities [m/s]

n_steps = numel(traj.t);
dt = traj.dt;
m  = size(rp0, 2);

r = zeros(3, m, n_steps);
r(:, :, 1) = rp0;
rk = rp0;
vk = vp0;
for k = 1:n_steps-1
    GM      = traj.GM(k, :);
    rb_k    = traj.r(:, :, k);
    rb_half = 0.5 * (traj.r(:, :, k) + traj.r(:, :, k+1));
    a1      = test_particle_acceleration(rk, rb_k, GM);
    r_half  = rk + vk * (dt / 2);
    v_half  = vk + a1 * (dt / 2);
    a2      = test_particle_acceleration(r_half, rb_half, GM);
    rk      = rk + v_half * dt;
    vk      = vk + a2 * dt;
    r(:, :, k+1) = rk;
end

part.t = traj.t;
part.r = r;
part.v = vk;
end
