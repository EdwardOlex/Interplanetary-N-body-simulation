function run_tests()
%RUN_TESTS  Verification and validation suite. Runs in MATLAB and GNU Octave.
%
%   From the repository root:   addpath tests; run_tests

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

tests = {@test_acceleration_matches_pairwise_sum, ...
         @test_circular_orbit_closes, ...
         @test_conservation_laws, ...
         @test_earth_against_jpl_horizons, ...
         @test_jupiter_against_jpl_horizons, ...
         @test_adaptive_integrator, ...
         @test_substep_impact_detection, ...
         @test_identical_trajectories_do_not_deviate};

n_fail = 0;
for i = 1:numel(tests)
    name = func2str(tests{i});
    try
        tic; msg = tests{i}(); el = toc;
        fprintf('  PASS  %-45s %6.1f s  %s\n', name, el, msg);
    catch err
        n_fail = n_fail + 1;
        fprintf('  FAIL  %-45s %s\n', name, err.message);
    end
end
fprintf('\n%d of %d tests passed\n', numel(tests) - n_fail, numel(tests));
if n_fail > 0
    error('run_tests:failed', '%d test(s) failed', n_fail);
end
end


function check(cond, varargin)
if ~cond, error(varargin{:}); end
end


function msg = test_acceleration_matches_pairwise_sum()
rng(1);
r = randn(3, 7) * 1e11;
GM = rand(1, 7) * 1e20;
a = nbody_acceleration(r, GM);
a_ref = zeros(3, 7);
for i = 1:7
    for j = 1:7
        if i ~= j
            d = r(:, i) - r(:, j);
            a_ref(:, i) = a_ref(:, i) - GM(j) * d / norm(d)^3;
        end
    end
end
err = max(abs(a(:) - a_ref(:))) / max(abs(a_ref(:)));
check(err < 1e-12, 'relative error %.2e', err);
msg = sprintf('rel. error %.1e', err);
end


function msg = test_circular_orbit_closes()
c = constants();
sys.names = {'Star', 'Planet'};
sys.GM = [c.G * 1.98847e30, 0];
sys.mass = sys.GM / c.G;
sys.radius = [0 0];
sys.r0 = [0 0; 0 0; 0 0]; sys.r0(1, 2) = c.AU;
vc = sqrt(sys.GM(1) / c.AU);
sys.v0 = [0 0; 0 vc; 0 0];
sys.epoch = 0;
P = 2 * pi * sqrt(c.AU^3 / sys.GM(1));
traj = propagate_nbody(sys, P, 20001);
err = norm(traj.r(:, 2, end) - sys.r0(:, 2)) / c.AU;
check(err < 1e-5, 'orbit does not close: %.2e au', err);
msg = sprintf('closure error %.1e au', err);
end


function msg = test_conservation_laws()
c = constants();
sys = solar_system_initial_state();
traj = propagate_nbody(sys, 2 * c.YEAR, round(2 * 365.25 * 24) + 1);
d = conservation_diagnostics(traj);
dp = max(sqrt(sum((d.momentum - d.momentum(:, 1)).^2, 1))) / max(sqrt(sum(d.momentum.^2, 1)) + eps);
check(d.energy_drift < 1e-6, 'energy drift %.2e', d.energy_drift);
check(d.L_drift < 1e-9, 'angular momentum drift %.2e', d.L_drift);
msg = sprintf('dE/E %.1e, dL/L %.1e, dP %.1e', d.energy_drift, d.L_drift, dp);
end


function msg = test_earth_against_jpl_horizons()
% Earth (399) barycentric ecliptic J2000 position on 2020-01-01 00:00 TDB, JPL Horizons
c = constants();
sys = solar_system_initial_state();
traj = propagate_nbody(sys, 365 * c.DAY, 365 * 24 + 1);
ref = [-1.701450916047793e-01; 9.765600183819082e-01; -1.813218493670265e-05] * c.AU;
err = norm(traj.r(:, 4, end) - ref);
check(err < 1000e3, 'position error %.0f km', err / 1e3);
msg = sprintf('1-year error %.0f km', err / 1e3);
end


function msg = test_jupiter_against_jpl_horizons()
% Jupiter (599) barycentric ecliptic J2000 position on 2020-01-01 00:00 TDB, JPL Horizons
c = constants();
sys = solar_system_initial_state();
traj = propagate_nbody(sys, 365 * c.DAY, 365 * 4 + 1);
ref = [5.223477248448145e-01; -5.193582882397713e+00; 9.853626938394582e-03] * c.AU;
err = norm(traj.r(:, 6, end) - ref);
check(err < 1e5 * 1e3, 'position error %.0f km', err / 1e3);
msg = sprintf('1-year error %.0f km (%.1e of its distance)', err / 1e3, err / norm(ref));
end


function msg = test_adaptive_integrator()
c = constants();
sys = solar_system_initial_state();
t_end = 30 * c.DAY;
a = propagate_nbody(sys, t_end, 31, [], 'ode45');
b = propagate_nbody(sys, t_end, 30 * 96 + 1);           % 15-minute fixed step
diff = max(sqrt(sum((a.r(:, :, end) - b.r(:, :, end)).^2, 1)));
check(diff < 50e3, 'ode45 and midpoint differ by %.0f km', diff / 1e3);
msg = sprintf('max difference %.2f km after 30 days', diff / 1e3);
end


function msg = test_substep_impact_detection()
% A fast particle crosses a planet between two samples of the time grid.
traj.t = [0 3600 7200];
traj.dt = 3600;
traj.names = {'Planet'};
traj.radius = 6.371e6;
traj.r = zeros(3, 1, 3);
traj.GM = zeros(3, 1);
v = 1e5;                                  % 100 km/s, i.e. 360 000 km per step
part.t = traj.t;
part.r = zeros(3, 1, 3);
part.r(:, 1, 1) = [-1.8e8; 1e6; 0];
part.r(:, 1, 2) = [-1.8e8 + v * 3600; 1e6; 0];
part.r(:, 1, 3) = [-1.8e8 + 2 * v * 3600; 1e6; 0];
tg.name = 'Planet'; tg.zones = 1e7; tg.labels = {'warning'};
enc = close_approach_screening(traj, part, tg);
check(enc.n_entries(end) == 1, 'impact not detected');
check(abs(enc.miss_distance - 1e6) < 1, 'miss distance %.1f m', enc.miss_distance);
t_entry = (1.8e8 - sqrt(6.371e6^2 - 1e6^2)) / v;
check(abs(enc.entry_time(end) - t_entry) < 1e-3, 'entry time off by %.3g s', enc.entry_time(end) - t_entry);
msg = 'impact found between samples';
end


function msg = test_identical_trajectories_do_not_deviate()
c = constants();
sys = solar_system_initial_state();
a = propagate_nbody(sys, 10 * c.DAY, 41);
cmp = compare_trajectories(a, a, 1);
check(all(isnan(cmp.first_exceed)) && all(cmp.max_deviation == 0), 'spurious deviation');
msg = 'no deviation';
end
