%EX02_STELLAR_FLYBY  A solar-mass intruder crossing the outer Solar System.
%
%   Adds a body of one solar mass at 20 au, moving at (-4e-3, 4e-3, 0)
%   au/day, and follows how the planetary orbits respond over ten years.
%   The specific angular momentum of each planet about the Sun is plotted:
%   it stays constant for an unperturbed Keplerian orbit (Kepler's second
%   law), so any change measures the strength of the perturbation.
%
%   Jupiter passes within about 0.03 au of the intruder, an encounter far
%   too close for a fixed step of hours. This example therefore uses the
%   adaptive integrator ('ode45'), which keeps the relative energy error
%   below 1e-10. It takes a few minutes.

clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
c = constants();

years = 10;
output_step_days = 1;

sys = solar_system_initial_state();
sys = add_body(sys, 'Intruder', sys.mass(1), [20 0 0], [-4e-3 4e-3 0]);

n_steps = round(years * 365.25 / output_step_days) + 1;
traj = propagate_nbody(sys, years * c.YEAR, n_steps, [], 'ode45');

plot_orbits(traj);
title('Solar-mass flyby');

d = conservation_diagnostics(traj);
figure('Color', 'w'); hold on
planets = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
for i = 1:numel(planets)
    j = find(strcmp(traj.names, planets{i}));
    col = body_style(planets{i});
    plot(traj.t / c.YEAR, d.h_sun(j, :) / d.h_sun(j, 1), 'Color', col, 'LineWidth', 1.2);
end
xlabel('time [yr]'); ylabel('h(t) / h(0)');
legend(planets, 'Location', 'best'); grid on; box on
title('Specific angular momentum about the Sun');
fprintf('Relative energy drift: %.2e\n', d.energy_drift);
