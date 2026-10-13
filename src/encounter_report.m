function tbl = encounter_report(enc, epoch, n_particles, list_events)
%ENCOUNTER_REPORT  Print a close-approach and impact summary.
%
%   ENCOUNTER_REPORT(enc, epoch, n_particles) prints, for every target in
%   ENC (CLOSE_APPROACH_SCREENING), how many of the N_PARTICLES simulated
%   asteroids entered each alert zone, the impact probability with its
%   binomial 1-sigma uncertainty (or the 95% upper limit 3/N when no
%   impact occurs) and the epoch of the first impact.
%
%   ENCOUNTER_REPORT(..., true) also lists every event with its date.
%
%   tbl = ENCOUNTER_REPORT(...) returns the summary as a struct array.

if nargin < 4, list_events = false; end
c = constants();
tbl = struct([]);

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf(' CLOSE-APPROACH SCREENING  (%d asteroids, epoch %s TDB)\n', ...
        n_particles, datestr(epoch, 'yyyy-mm-dd HH:MM'));
fprintf('%s\n', repmat('=', 1, 78));

for i = 1:numel(enc)
    e = enc(i);
    nz = numel(e.zones);
    n_imp = e.n_entries(end);
    p  = n_imp / n_particles;
    sp = sqrt(max(p * (1 - p), 0) / n_particles);
    fprintf('\n %s\n', upper(e.name));
    for z = 1:nz
        if strcmp(e.labels{z}, 'impact')
            fprintf('   %-9s (r < %9.0f km)   : %5d\n', e.labels{z}, e.zones(z) / 1e3, e.n_entries(z));
        else
            fprintf('   %-9s (d < %9.3f au)   : %5d\n', e.labels{z}, e.zones(z) / c.AU, e.n_entries(z));
        end
    end
    if n_imp > 0
        fprintf('   impact probability             : %.4f +/- %.4f\n', p, sp);
    else
        fprintf('   impact probability             : < %.4f (95%% upper limit, no impacts)\n', 3 / n_particles);
    end
    [dmin, imin] = min(e.miss_distance);
    fprintf('   closest approach               : %.4g km (asteroid %d, %s)\n', ...
            dmin / 1e3, imin, epoch_to_date(epoch, e.tca(imin)));
    first_imp = min(e.entry_time(end, :));
    if ~isnan(first_imp)
        fprintf('   first impact                   : %s\n', epoch_to_date(epoch, first_imp));
    end

    if list_events
        for z = 1:nz
            ids = find(~isnan(e.entry_time(z, :)));
            if isempty(ids), continue; end
            [tt, order] = sort(e.entry_time(z, ids));
            ids = ids(order);
            fprintf('     %s events:\n', e.labels{z});
            for k = 1:numel(ids)
                fprintf('       asteroid %5d   %s\n', ids(k), epoch_to_date(epoch, tt(k)));
            end
        end
    end

    tbl(i).name = e.name;
    tbl(i).labels = e.labels;
    tbl(i).n_entries = e.n_entries;
    tbl(i).impact_probability = p;
    tbl(i).impact_probability_sigma = sp;
end
fprintf('\n');
end
