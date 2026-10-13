function cmp = compare_trajectories(ref, pert, threshold)
%COMPARE_TRAJECTORIES  Detect orbital deviations between a nominal and a perturbed propagation.
%
%   cmp = COMPARE_TRAJECTORIES(ref, pert, threshold) compares two
%   trajectories produced by PROPAGATE_NBODY on the same time grid, for
%   example the nominal Solar System and the same system with an extra
%   perturbing body. Only the bodies present in REF are compared.
%
%   For each body it returns the position deviation |r_pert - r_ref| along
%   the propagation and the first epoch at which it exceeds THRESHOLD [m].
%
%   Output struct CMP
%     names           1xN body names
%     deviation       NxT position deviation [m]
%     max_deviation   1xN maximum deviation [m]
%     first_exceed    1xN first epoch with deviation > threshold [s], NaN if never
%     threshold       threshold used [m]

if numel(ref.t) ~= numel(pert.t) || any(abs(ref.t - pert.t) > 1e-6 * ref.dt)
    error('compare_trajectories:grid', 'Both trajectories must share the same time grid.');
end
nb = size(ref.r, 2);
dev = reshape(sqrt(sum((pert.r(:, 1:nb, :) - ref.r).^2, 1)), nb, []);
first = nan(1, nb);
for i = 1:nb
    k = find(dev(i, :) > threshold, 1);
    if ~isempty(k), first(i) = ref.t(k); end
end

cmp.names = ref.names;
cmp.deviation = dev;
cmp.max_deviation = max(dev, [], 2).';
cmp.first_exceed = first;
cmp.threshold = threshold;
end
