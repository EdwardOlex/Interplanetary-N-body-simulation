function ax = plot_orbits(traj, bodies, ax)
%PLOT_ORBITS  3D plot of propagated orbits, with the final position of each body.
%
%   PLOT_ORBITS(traj) plots every body of TRAJ (PROPAGATE_NBODY), in au.
%   PLOT_ORBITS(traj, bodies) plots only the bodies named in the cell array BODIES.
%   PLOT_ORBITS(traj, bodies, ax) draws on the axes AX.

c = constants();
if nargin < 2 || isempty(bodies), bodies = traj.names; end
if nargin < 3 || isempty(ax), figure('Color', 'w'); ax = gca; end
hold(ax, 'on');
for i = 1:numel(bodies)
    j = find(strcmp(traj.names, bodies{i}), 1);
    if isempty(j), continue; end
    [col, ms] = body_style(bodies{i});
    xyz = reshape(traj.r(:, j, :), 3, []) / c.AU;
    plot3(ax, xyz(1, :), xyz(2, :), xyz(3, :), '-', 'Color', col, 'LineWidth', 1);
    plot3(ax, xyz(1, end), xyz(2, end), xyz(3, end), '.', 'Color', col, 'MarkerSize', ms);
end
xlabel(ax, 'x [au]'); ylabel(ax, 'y [au]'); zlabel(ax, 'z [au]');
axis(ax, 'equal'); grid(ax, 'on'); box(ax, 'on');
view(ax, 30, 25);
hold(ax, 'off');
end
