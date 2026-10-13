function enc = close_approach_screening(traj, part, targets)
%CLOSE_APPROACH_SCREENING  Close-approach and impact screening of test particles against target bodies.
%
%   enc = CLOSE_APPROACH_SCREENING(traj, part, targets) compares the
%   trajectories of the test particles PART (PROPAGATE_TEST_PARTICLES) with
%   those of the massive bodies in TRAJ (PROPAGATE_NBODY).
%
%   TARGETS is a struct array with fields
%     name     body name, as in traj.names
%     zones    1xZ radii of the alert zones around the body [m]
%     labels   1xZ cell array with a label for each zone
%   An impact zone equal to the physical radius of the body is always
%   added as the last zone, with label 'impact'.
%
%   Between two ephemeris epochs the relative motion is assumed linear, so
%   encounters are resolved below the integration step: a fast asteroid
%   that crosses a planet between two samples is still detected.
%
%   Output struct array ENC (one element per target)
%     name           target name
%     miss_distance  1xM minimum distance to the target [m]
%     tca            1xM time of closest approach [s since epoch]
%     zones, labels  alert zones (impact included)
%     entry_time     ZxM time of first entry into each zone [s], NaN if never
%     n_entries      1xZ number of particles that entered each zone

T  = numel(traj.t);
dt = traj.dt;
m  = size(part.r, 2);
enc = struct([]);

for it = 1:numel(targets)
    tg = targets(it);
    j = find(strcmp(traj.names, tg.name), 1);
    if isempty(j)
        error('close_approach_screening:unknownTarget', 'Unknown target "%s".', tg.name);
    end
    zones  = [tg.zones(:).', traj.radius(j)];
    labels = [tg.labels(:).', {'impact'}];

    nz = numel(zones);
    miss  = zeros(m, 1);
    tca   = zeros(m, 1);
    entry = nan(nz, m);

    % Process the particles in chunks to bound memory use
    chunk = 100;
    for i0 = 1:chunk:m
        ids = i0:min(i0 + chunk - 1, m);
        mc  = numel(ids);
        rel = part.r(:, ids, :) - traj.r(:, j, :);     % 3 x mc x T
        d0  = rel(:, :, 1:T-1);
        dd  = rel(:, :, 2:T) - d0;
        A = reshape(sum(dd.^2, 1),    mc, T-1);
        B = reshape(sum(d0 .* dd, 1), mc, T-1);
        C = reshape(sum(d0.^2, 1),    mc, T-1);
        clear rel d0 dd

        % Minimum distance within each step, relative motion linear in time
        s = min(max(-B ./ A, 0), 1);
        s(A == 0) = 0;
        dmin = sqrt(max(C + 2 * B .* s + A .* s.^2, 0));   % mc x (T-1)

        [mc_miss, kmin] = min(dmin, [], 2);
        lin = sub2ind([mc, T-1], (1:mc).', kmin);
        miss(ids) = mc_miss;
        tca(ids)  = traj.t(kmin).' + s(lin) * dt;

        % First entry into each zone: smallest root of |d0 + s dd| = R
        for z = 1:nz
            Rz = zones(z);
            [hit, kf] = max(dmin < Rz, [], 2);
            loc = find(hit);
            if isempty(loc), continue; end
            lin = sub2ind([mc, T-1], loc, kf(loc));
            a = A(lin); b = B(lin); cc = C(lin);
            disc = max(b.^2 - a .* (cc - Rz^2), 0);
            se = (-b - sqrt(disc)) ./ a;
            se(cc <= Rz^2 | a == 0) = 0;
            se = min(max(se, 0), 1);
            entry(z, ids(loc)) = traj.t(kf(loc)) + (se.') * dt;
        end
    end

    enc(it).name          = tg.name;
    enc(it).miss_distance = miss.';
    enc(it).tca           = tca.';
    enc(it).zones         = zones;
    enc(it).labels        = labels;
    enc(it).entry_time    = entry;
    enc(it).n_entries     = sum(~isnan(entry), 2).';
end
end
