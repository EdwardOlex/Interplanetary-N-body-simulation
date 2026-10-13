function ax = plot_impact_assessment(traj, part, enc, target, ax)
%PLOT_IMPACT_ASSESSMENT  Asteroid trajectories coloured by their encounter with a target body.
%
%   PLOT_IMPACT_ASSESSMENT(traj, part, enc, target) plots the planetary
%   orbits and the asteroid tracks, highlighting the asteroids that entered
%   the alert zones of TARGET (name of a body in ENC, default 'Earth'):
%   grey: no alert, orange: warning, red: critical or impact.

c = constants();
if nargin < 4 || isempty(target), target = 'Earth'; end
if nargin < 5 || isempty(ax), figure('Color', 'w'); ax = gca; end
e = enc(strcmp({enc.name}, target));
if isempty(e), error('plot_impact_assessment:target', 'Target "%s" not screened.', target); end

level = zeros(1, size(part.r, 2));
for z = 1:numel(e.zones)
    level(~isnan(e.entry_time(z, :))) = z;
end
cols = [0.70 0.70 0.70; 0.95 0.55 0.10; 0.85 0.10 0.10; 0.55 0.00 0.00];

hold(ax, 'on');
for lv = 0:numel(e.zones)
    ids = find(level == lv);
    for i = ids
        xyz = reshape(part.r(:, i, :), 3, []) / c.AU;
        lw = 0.5 + 0.6 * (lv > 0);
        plot3(ax, xyz(1, :), xyz(2, :), xyz(3, :), '-', 'Color', cols(min(lv, 3) + 1, :), 'LineWidth', lw);
    end
end
plot_orbits(traj, {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter'}, ax);
title(ax, sprintf('%s: %d warning, %d critical, %d impact', target, e.n_entries(1), ...
      e.n_entries(min(2, end)), e.n_entries(end)));
hold(ax, 'off');
end
