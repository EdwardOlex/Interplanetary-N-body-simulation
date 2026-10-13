function [rp0, vp0] = sample_asteroid_cloud(n, shell_radius, speed, max_angle)
%SAMPLE_ASTEROID_CLOUD  Monte Carlo initial conditions for an inbound population of asteroids.
%
%   [rp0, vp0] = SAMPLE_ASTEROID_CLOUD(n, shell_radius, speed) draws N
%   asteroids uniformly distributed on a sphere of radius SHELL_RADIUS [m]
%   centred on the origin, all moving inwards with speed SPEED [m/s].
%
%   The velocity of each asteroid starts pointing to the centre of the
%   sphere and is then rotated by three random angles about the x, y and z
%   axes, each drawn uniformly in [-max_angle, max_angle] (default pi/4).
%   This produces a broad distribution of inbound trajectories, most of
%   which cross the inner Solar System.
%
%   Set the random seed beforehand (rng(seed) in MATLAB) for reproducible
%   samples.

if nargin < 4, max_angle = pi / 4; end

% Uniform on the sphere: azimuth uniform in [0, 2 pi), cos(polar) uniform in [-1, 1]
phi = 2 * pi * rand(1, n);
cos_theta = -1 + 2 * rand(1, n);
sin_theta = sqrt(1 - cos_theta.^2);
rp0 = shell_radius * [sin_theta .* cos(phi); sin_theta .* sin(phi); cos_theta];

% Inbound direction, randomly tilted
inward = -rp0 ./ sqrt(sum(rp0.^2, 1));
vp0 = zeros(3, n);
for i = 1:n
    ax = -max_angle + 2 * max_angle * rand;
    ay = -max_angle + 2 * max_angle * rand;
    az = -max_angle + 2 * max_angle * rand;
    Rx = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)];
    Ry = [cos(ay) 0 sin(ay); 0 1 0; -sin(ay) 0 cos(ay)];
    Rz = [cos(az) -sin(az) 0; sin(az) cos(az) 0; 0 0 1];
    vp0(:, i) = speed * (Rx * Ry * Rz * inward(:, i));
end
end
