function animate_orbits(traj, filename, part, lim_au, every)
%ANIMATE_ORBITS  Write a video of the propagated system.
%
%   ANIMATE_ORBITS(traj, filename) writes the motion of all bodies of TRAJ
%   to the video FILENAME (MATLAB VideoWriter).
%   ANIMATE_ORBITS(traj, filename, part) also shows the test particles PART.
%   ANIMATE_ORBITS(..., lim_au, every) sets the half-width of the box [au]
%   (default 30) and writes one frame every EVERY epochs (default 10).

c = constants();
if nargin < 3, part = []; end
if nargin < 4 || isempty(lim_au), lim_au = 30; end
if nargin < 5 || isempty(every), every = 10; end
if ~exist('VideoWriter', 'class') && ~exist('VideoWriter', 'file')
    warning('animate_orbits:noVideo', 'VideoWriter is not available; no video written.');
    return
end

fig = figure('Color', 'w');
vid = VideoWriter(filename);
open(vid);
nb = size(traj.r, 2);
for k = 1:every:numel(traj.t)
    clf(fig); hold on
    for j = 1:nb
        [col, ms] = body_style(traj.names{j});
        p = traj.r(:, j, k) / c.AU;
        plot3(p(1), p(2), p(3), '.', 'Color', col, 'MarkerSize', ms);
    end
    if ~isempty(part)
        p = part.r(:, :, k) / c.AU;
        plot3(p(1, :), p(2, :), p(3, :), '.k', 'MarkerSize', 4);
    end
    axis equal; xlim(lim_au * [-1 1]); ylim(lim_au * [-1 1]); zlim(lim_au * [-1 1]);
    xlabel('x [au]'); ylabel('y [au]'); zlabel('z [au]'); view(30, 25); grid on
    title(epoch_to_date(traj.epoch, traj.t(k)));
    hold off
    writeVideo(vid, getframe(fig));
end
close(vid);
close(fig);
end
