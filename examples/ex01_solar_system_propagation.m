%EX01_SOLAR_SYSTEM_PROPAGATION  Nominal propagation of the Sun, the eight planets and the Moon.
%
%   Propagates the Solar System for ten years from 2019-01-01 TDB, plots
%   the orbits and reports the conservation of energy and angular
%   momentum. Optionally writes an animation (MATLAB only).

clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
c = constants();

% ---- configuration ---------------------------------------------------
years      = 10;          % propagation span [yr]
step_hours = 6;           % integration step [h]
make_video = false;       % write SolarSystem.avi (MATLAB only)

% ---- propagation -----------------------------------------------------
sys = solar_system_initial_state();
n_steps = round(years * 365.25 * 24 / step_hours) + 1;
tic
traj = propagate_nbody(sys, years * c.YEAR, n_steps);
fprintf('Propagated %d bodies over %g years (%d steps of %.2f h) in %.1f s\n', ...
        numel(sys.names), years, n_steps - 1, traj.dt / 3600, toc);

% ---- validation ------------------------------------------------------
d = conservation_diagnostics(traj);
fprintf('Relative energy drift            : %.2e\n', d.energy_drift);
fprintf('Relative angular-momentum drift  : %.2e\n', d.L_drift);

% ---- figures ---------------------------------------------------------
plot_orbits(traj);
title(sprintf('Solar System, %s to %s TDB', epoch_to_date(traj.epoch, 0), ...
      epoch_to_date(traj.epoch, traj.t(end))));

figure('Color', 'w');
plot_orbits(traj, {'Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars'}, gca);
title('Inner Solar System');

if make_video
    animate_orbits(traj, 'SolarSystem', [], 32, 20);
end
