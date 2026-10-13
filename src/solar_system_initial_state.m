function sys = solar_system_initial_state()
%SOLAR_SYSTEM_INITIAL_STATE  Initial state vectors of the Sun, planets and Moon.
%
%   sys = SOLAR_SYSTEM_INITIAL_STATE() returns a struct with the barycentric
%   state of the Solar System at the epoch 2019-01-01 00:00:00 TDB.
%
%   Source of the state vectors: JPL Horizons (https://ssd.jpl.nasa.gov/horizons/),
%   centre Solar System Barycentre, reference frame ecliptic of J2000.0.
%
%   Fields
%     names   1xN cell array of body names
%     GM      1xN gravitational parameters [m^3 s^-2]
%     mass    1xN masses [kg] (GM / G)
%     radius  1xN mean radii [m]
%     r0      3xN positions [m]
%     v0      3xN velocities [m s^-1]
%     epoch   epoch as a datenum (TDB)
%
%   Additional bodies can be appended with ADD_BODY.

c = constants();

names = {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', ...
         'Saturn', 'Uranus', 'Neptune', 'Moon'};

% Positions [au] and velocities [au/day], 2019-01-01 00:00 TDB
r_au = [ -7.558735329017875e-04   7.428856011015974e-03  -5.737173235978428e-05   % Sun
         -3.159243724530365e-01  -3.122664112956987e-01   2.732130737118980e-03   % Mercury
         -5.506840617208332e-01   4.682823885525094e-01   3.800118466168461e-02   % Venus
         -1.715115335386781e-01   9.758005151474849e-01  -1.036379644680106e-04   % Earth
          1.089458162138731e+00   9.687105351565151e-01  -6.666116417766668e-03   % Mars
         -2.135047980135147e+00  -4.897850912579155e+00   6.807194735366415e-02   % Jupiter
          1.958822848902041e+00  -9.859694810154547e+00   9.346113213197697e-02   % Saturn
          1.701531007391380e+01   1.025052953139990e+01  -1.823648889669433e-01   % Uranus
          2.898163620691469e+01  -7.481208033523555e+00  -5.138498122721312e-01   % Neptune
         -1.734238606952717e-01   9.740726777725756e-01   1.255227657142639e-04]; % Moon

v_au = [ -7.859449720803353e-06   1.979914352132652e-06   2.008516510572875e-07
          1.430292509441845e-02  -1.846637695271600e-02  -2.821753137305481e-03
         -1.307514440770031e-02  -1.560378950041907e-02   5.401533925256370e-04
         -1.723595211516194e-02  -3.045322548971216e-03   5.346408248909461e-07
         -8.731137654316889e-03   1.169226996727862e-02   4.592058136523131e-04
          6.826434975087506e-03  -2.654951723393821e-03  -1.416782069788317e-04
          5.164140819216792e-03   1.070076264025756e-03  -2.244505966450992e-04
         -2.058436023906124e-03   3.185655137014581e-03   3.844838369437882e-05
          7.638654639184098e-04   3.057829840918665e-03  -8.097651119274648e-05
         -1.686465147855430e-02  -3.498039880324010e-03  -1.136004747701177e-05];

% Gravitational parameters [km^3 s^-2] (JPL DE440; planets include their
% satellites except for the Earth and the Moon, which are separate bodies)
GM_km3 = [1.32712440041279e11, 2.2031868551e4, 3.24858592e5, 3.98600435507e5, ...
          4.2828375816e4, 1.267127641e8, 3.79405848418e7, 5.7945564e6, ...
          6.8365271005e6, 4.902800118e3];
% Mean radii [m]
radius = [6.957e8, 2.4397e6, 6.0518e6, 6.371e6, 3.3895e6, ...
          6.9911e7, 5.8232e7, 2.5362e7, 2.4622e7, 1.7374e6];

sys.names  = names;
sys.GM     = GM_km3 * 1e9;
sys.mass   = sys.GM / c.G;
sys.radius = radius;
sys.r0     = (r_au * c.AU).';
sys.v0     = (v_au * c.AU / c.DAY).';
sys.epoch  = datenum(2019, 1, 1, 0, 0, 0);
end
