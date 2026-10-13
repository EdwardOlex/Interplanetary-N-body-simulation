%EX05_ASTEROID_IMPACT_MONTE_CARLO  Monte Carlo close-approach and impact assessment for an inbound asteroid population.
%
%   1. The Sun, planets and Moon are propagated from 2019-01-01 TDB.
%   2. A population of asteroids is drawn on a sphere of 5 au and sent
%      inwards with random inclinations (Monte Carlo sampling).
%   3. The asteroids are propagated as test particles through the
%      planetary ephemeris (restricted N-body problem).
%   4. Every asteroid is screened against the planets: alert zones,
%      miss distance, time of closest approach and impact, with sub-step
%      resolution of each encounter.
%   5. Impact probabilities and their binomial uncertainties are reported.

clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
c = constants();

% ---- configuration ---------------------------------------------------
n_asteroids  = 1000;                  % Monte Carlo sample size
shell_radius = 5 * c.AU;              % initial distance from the barycentre
speed        = 3e-3 * c.AU / c.DAY;   % initial speed (about 5.2 km/s)
years        = 7;                     % propagation span
step_hours   = 6;                     % integration step
seed         = 2019;                  % random seed
list_events  = false;                 % print every event with its date

targets = planetary_targets({'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter'});

% ---- propagation -----------------------------------------------------
rng(seed);
sys = solar_system_initial_state();
n_steps = round(years * 365.25 * 24 / step_hours) + 1;

tic
traj = propagate_nbody(sys, years * c.YEAR, n_steps);
[rp0, vp0] = sample_asteroid_cloud(n_asteroids, shell_radius, speed);
part = propagate_test_particles(traj, rp0, vp0);
fprintf('Propagated %d asteroids over %g years in %.1f s\n', n_asteroids, years, toc);

% ---- screening -------------------------------------------------------
enc = close_approach_screening(traj, part, targets);
encounter_report(enc, traj.epoch, n_asteroids, list_events);

% ---- figures ---------------------------------------------------------
plot_impact_assessment(traj, part, enc, 'Earth');

figure('Color', 'w');
e = enc(strcmp({enc.name}, 'Earth'));
histogram_au = e.miss_distance / c.AU;
edges = logspace(log10(max(min(histogram_au), 1e-5)), log10(max(histogram_au)), 40);
counts = histc(histogram_au, edges);
stairs(edges, counts, 'LineWidth', 1.2); set(gca, 'XScale', 'log'); grid on
xlabel('miss distance to the Earth [au]'); ylabel('asteroids');
title('Distribution of Earth miss distances');
