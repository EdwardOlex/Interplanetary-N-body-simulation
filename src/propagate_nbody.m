function traj = propagate_nbody(sys, t_end, n_steps, gm_schedule, method)
%PROPAGATE_NBODY  Propagate a gravitating N-body system.
%
%   traj = PROPAGATE_NBODY(sys, t_end, n_steps) integrates the system SYS
%   (see SOLAR_SYSTEM_INITIAL_STATE) from t = 0 to t = T_END [s] on a fixed
%   grid of N_STEPS epochs (time step dt = t_end / (n_steps - 1)).
%
%   traj = PROPAGATE_NBODY(..., gm_schedule) lets the gravitational
%   parameters change during the propagation. GM_SCHEDULE is a function
%   handle @(k, t) returning a 1xN vector of GM [m^3 s^-2] to be used
%   during step k (t is the epoch at the start of the step). This is used,
%   for example, to model a growing Jupiter.
%
%   traj = PROPAGATE_NBODY(..., method) selects the integrator:
%     'midpoint' (default) fixed-step, second-order explicit midpoint scheme
%                  r_half = r_k + v_k dt/2,   v_half = v_k + a(r_k) dt/2
%                  r_k+1  = r_k + v_half dt,  v_k+1  = v_k + a(r_half) dt
%     'ode45'    adaptive Dormand-Prince 4(5) with relative tolerance 1e-11,
%                for close encounters that a fixed step cannot resolve.
%                The state is still returned on the fixed output grid.
%   With 'ode45', GM_SCHEDULE is called as gm_schedule(NaN, t).
%
%   Output struct TRAJ
%     t      1xT epochs since SYS.epoch [s]
%     r, v   3xNxT positions [m] and velocities [m/s]
%     GM     TxN gravitational parameters used at each epoch
%     names, radius, epoch   copied from SYS

if nargin < 4, gm_schedule = []; end
if nargin < 5 || isempty(method), method = 'midpoint'; end

nb = size(sys.r0, 2);
dt = t_end / (n_steps - 1);
t  = (0:n_steps-1) * dt;

r = zeros(3, nb, n_steps);
v = zeros(3, nb, n_steps);
GM_hist = repmat(sys.GM(:).', n_steps, 1);
r(:, :, 1) = sys.r0;
v(:, :, 1) = sys.v0;

if strcmpi(method, 'ode45')
    [r, v, GM_hist] = propagate_ode45(sys, t, gm_schedule);
    traj = pack_traj(sys, t, dt, r, v, GM_hist);
    return
elseif ~strcmpi(method, 'midpoint')
    error('propagate_nbody:method', 'Unknown method "%s".', method);
end

rk = sys.r0;
vk = sys.v0;
GM = sys.GM;
for k = 1:n_steps-1
    if ~isempty(gm_schedule)
        GM = gm_schedule(k, t(k));
        GM_hist(k, :) = GM;
    end
    a1     = nbody_acceleration(rk, GM);
    r_half = rk + vk * (dt / 2);
    v_half = vk + a1 * (dt / 2);
    a2     = nbody_acceleration(r_half, GM);
    rk     = rk + v_half * dt;
    vk     = vk + a2 * dt;
    r(:, :, k+1) = rk;
    v(:, :, k+1) = vk;
end
if ~isempty(gm_schedule)
    GM_hist(end, :) = GM;
end

traj = pack_traj(sys, t, dt, r, v, GM_hist);
end


function traj = pack_traj(sys, t, dt, r, v, GM_hist)
traj.t      = t;
traj.dt     = dt;
traj.r      = r;
traj.v      = v;
traj.GM     = GM_hist;
traj.names  = sys.names;
traj.radius = sys.radius;
traj.epoch  = sys.epoch;
end


function [r, v, GM_hist] = propagate_ode45(sys, t, gm_schedule)
nb = size(sys.r0, 2);
if isempty(gm_schedule)
    gmf = @(tt) sys.GM;
else
    gmf = @(tt) gm_schedule(NaN, tt);
end
rhs = @(tt, y) [y(3*nb+1:end); reshape(nbody_acceleration(reshape(y(1:3*nb), 3, nb), gmf(tt)), [], 1)];
y0 = [sys.r0(:); sys.v0(:)];
scale = max(abs(y0(1:3*nb)));
opts = odeset('RelTol', 1e-11, 'AbsTol', [1e-11 * scale * ones(3*nb, 1); 1e-11 * ones(3*nb, 1)]);
[~, Y] = ode45(rhs, t, y0, opts);
n = numel(t);
r = reshape(Y(:, 1:3*nb).', 3, nb, n);
v = reshape(Y(:, 3*nb+1:end).', 3, nb, n);
GM_hist = zeros(n, nb);
for k = 1:n
    GM_hist(k, :) = gmf(t(k));
end
end
