%EX03_PERTURBATION_SENSITIVITY  Sensitivity of the planetary ephemeris to a small perturbing body.
%
%   Propagates the nominal Solar System and the same system with an extra
%   body of one hundredth of Mercury's mass entering from 20 au. The two
%   ephemerides are compared epoch by epoch, and every body whose position
%   departs from the nominal one by more than a given tolerance is
%   reported with the date at which the deviation first exceeds it.

clear; close all; clc
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
c = constants();

years      = 10;
step_hours = 6;
tolerance  = 3e3;                      % position tolerance [m]

nominal = solar_system_initial_state();
perturbed = add_body(nominal, 'Perturber', nominal.mass(2) / 100, [20 0 0], [-4e-3 4e-3 0]);

n_steps = round(years * 365.25 * 24 / step_hours) + 1;
ref  = propagate_nbody(nominal,   years * c.YEAR, n_steps);
pert = propagate_nbody(perturbed, years * c.YEAR, n_steps);

cmp = compare_trajectories(ref, pert, tolerance);

fprintf('\nBodies displaced by more than %.0f km:\n', tolerance / 1e3);
for i = 1:numel(cmp.names)
    if ~isnan(cmp.first_exceed(i))
        fprintf('  %-8s  from %s   (max deviation %.1f km)\n', cmp.names{i}, ...
                epoch_to_date(ref.epoch, cmp.first_exceed(i)), cmp.max_deviation(i) / 1e3);
    end
end

figure('Color', 'w');
semilogy(ref.t / c.YEAR, max(cmp.deviation, 1e-3).'); grid on
xlabel('time [yr]'); ylabel('position deviation [m]');
legend(cmp.names, 'Location', 'southeast');
title('Deviation from the nominal ephemeris');
