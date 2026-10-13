function targets = planetary_targets(names)
%PLANETARY_TARGETS  Default alert zones for close-approach screening.
%
%   targets = PLANETARY_TARGETS(names) returns the target definition used by
%   CLOSE_APPROACH_SCREENING for the bodies listed in the cell array NAMES
%   (default {'Earth', 'Jupiter'}).
%
%   Each target has two alert zones besides the impact zone:
%     terrestrial planets and the Moon:  warning 0.30 au, critical 0.05 au
%     giant planets:                     warning 1.00 au, critical 0.30 au
%   The 0.05 au critical zone matches the minimum orbit intersection
%   distance used to classify potentially hazardous asteroids.

if nargin < 1, names = {'Earth', 'Jupiter'}; end
c = constants();
giants = {'Jupiter', 'Saturn', 'Uranus', 'Neptune'};

targets = struct('name', {}, 'zones', {}, 'labels', {});
for i = 1:numel(names)
    targets(i).name = names{i};
    if any(strcmp(names{i}, giants))
        targets(i).zones = [1.00, 0.30] * c.AU;
    else
        targets(i).zones = [0.30, 0.05] * c.AU;
    end
    targets(i).labels = {'warning', 'critical'};
end
end
