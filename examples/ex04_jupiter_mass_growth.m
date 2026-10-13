%EX04_JUPITER_MASS_GROWTH  Stability of the Solar System when Jupiter grows into a star.
%
%   Jupiter's mass is increased linearly from its present value to one
%   solar mass over the first 20 years and then kept constant. The
%   propagation spans 50 years and shows how the planetary system is
%   reshaped by a second star.
%
%   Once Jupiter becomes a star, the inner planets undergo repeated close
%   encounters and the evolution is chaotic: a fixed step of hours gives
%   results that change with the step size. The adaptive integrator
%   ('ode45') is therefore used. The run is slow, especially in GNU Octave.

clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
c = constants();

years       = 50;
output_step_days = 1;
ramp_years  = 20;

sys = solar_system_initial_state();
jup = find(strcmp(sys.names, 'Jupiter'));
GM0 = sys.GM;
GM_final = sys.GM(1);                         % one solar mass

n_steps = round(years * 365.25 / output_step_days) + 1;
only_jupiter = double(1:numel(GM0) == jup);
schedule = @(k, t) GM0 + only_jupiter * (GM_final - GM0(jup)) * min(t / (ramp_years * c.YEAR), 1);

traj = propagate_nbody(sys, years * c.YEAR, n_steps, schedule, 'ode45');

plot_orbits(traj);
title('Jupiter grows to one solar mass in 20 years');

