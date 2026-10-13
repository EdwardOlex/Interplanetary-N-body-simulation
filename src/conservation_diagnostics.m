function dg = conservation_diagnostics(traj)
%CONSERVATION_DIAGNOSTICS  Energy, momentum and angular momentum of a propagated N-body system.
%
%   dg = CONSERVATION_DIAGNOSTICS(traj) evaluates the integrals of motion
%   of the trajectory TRAJ (PROPAGATE_NBODY) at every epoch. For a closed
%   system with constant masses they are conserved by the exact dynamics,
%   so their drift measures the accuracy of the propagation.
%
%   Output struct DG
%     energy            1xT total energy [J]
%     angular_momentum  3xT total angular momentum [kg m^2 s^-1]
%     momentum          3xT total linear momentum [kg m s^-1]
%     energy_drift      max |E(t) - E(0)| / |E(0)|
%     L_drift           max |L(t) - L(0)| / |L(0)|
%     h_sun             NxT specific angular momentum of each body about the
%                       Sun [m^2 s^-1] (constant for a Keplerian orbit:
%                       Kepler's second law)

c = constants();
[~, nb, T] = size(traj.r);
E = zeros(1, T);
L = zeros(3, T);
P = zeros(3, T);
h = zeros(nb, T);
isun = find(strcmp(traj.names, 'Sun'), 1);

for k = 1:T
    m = traj.GM(k, :) / c.G;
    r = traj.r(:, :, k);
    v = traj.v(:, :, k);
    kin = 0.5 * sum(m .* sum(v.^2, 1));
    dx = r(1, :).' - r(1, :);
    dy = r(2, :).' - r(2, :);
    dz = r(3, :).' - r(3, :);
    d = sqrt(dx.^2 + dy.^2 + dz.^2);
    mm = m(:) * m(:).';
    pot = -c.G * sum(mm(triu(true(nb), 1)) ./ d(triu(true(nb), 1)));
    E(k) = kin + pot;
    L(:, k) = sum(m .* cross(r, v, 1), 2);
    P(:, k) = sum(m .* v, 2);
    if ~isempty(isun)
        rr = r - r(:, isun);
        vv = v - v(:, isun);
        h(:, k) = sqrt(sum(cross(rr, vv, 1).^2, 1)).';
    end
end

dg.energy = E;
dg.angular_momentum = L;
dg.momentum = P;
dg.energy_drift = max(abs(E - E(1))) / abs(E(1));
dg.L_drift = max(sqrt(sum((L - L(:, 1)).^2, 1))) / norm(L(:, 1));
dg.h_sun = h;
end
