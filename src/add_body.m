function sys = add_body(sys, name, mass, r_au, v_au_per_day, radius)
%ADD_BODY  Append a massive body to a system created by SOLAR_SYSTEM_INITIAL_STATE.
%
%   sys = ADD_BODY(sys, name, mass, r_au, v_au_per_day) adds a body of mass
%   MASS [kg] with initial position R_AU [au] and velocity V_AU_PER_DAY
%   [au/day], both 3-element vectors in the same frame as SYS.
%   sys = ADD_BODY(..., radius) also sets its radius [m] (default 0).
%
%   Example: a solar-mass intruder entering the outer Solar System
%     sys = add_body(sys, 'Intruder', 1.98847e30, [20 0 0], [-4e-3 4e-3 0]);

if nargin < 6, radius = 0; end
c = constants();
sys.names{end+1} = name;
sys.mass(end+1)  = mass;
sys.GM(end+1)    = c.G * mass;
sys.radius(end+1) = radius;
sys.r0(:, end+1) = r_au(:) * c.AU;
sys.v0(:, end+1) = v_au_per_day(:) * c.AU / c.DAY;
end
