function [x, y, z, vx, vy, vz] = orb2cart_ecc(a, ecc, inc, w, R, f)

%%  Coordinate transformation from orbital elements to cartesian coordinates.
%
%   Synopsis: this program transforms orbital elements (a, e, inc, w, R, f)
%           into cartesian elements in the ECI frame of reference.
%
%   INPUT
%       a: semimajor axis (km)
%       ecc: eccentricity
%       inc: inclination (deg)
%       w: argument of perigee (deg)
%       R: right ascension of the ascending node (deg)
%       f: mean anomaly (deg)
%
%   OUTPUT
%       x, y, z: cartesian position of the satellite (km)
%       vx, vy, vz: cartesian velocity of the satellite (km/s)
%
%   Authors and Updates
%       David Arnas 
%       26 January 2022

%% Earth gravitational constant
mu = 398600.4418;

%% Tranformation to radians
inc = inc*pi/180;
R = R*pi/180;
w = w*pi/180;
f = f*pi/180;

%% Transformation of coordinates (elliptic motion)  
% Semilatus rectum
p = a*(1 - ecc^2);

% Rotation matrix
g = [[cos(R)*cos(w) - sin(R)*sin(w)*cos(inc), ...
    -cos(R)*sin(w) - sin(R)*cos(w)*cos(inc), sin(R)*sin(inc)]; ...
    [sin(R)*cos(w) + cos(R)*sin(w)*cos(inc), ...
    -sin(R)*sin(w) + cos(R)*cos(w)*cos(inc), -cos(R)*sin(inc)]; ...
    [sin(w)*sin(inc), cos(w)*sin(inc), cos(inc)]];

% Coordinates in the perifocal frame of reference
xp = [cos(f), sin(f), 0]';
xp = p*xp/(1 + ecc*cos(f));
vp = [-sin(f), (ecc + cos(f)), 0]';
vp = sqrt(mu/p)*vp;

% Cartesian coordinates
xp = g*xp;
vp = g*vp;

x = xp(1);
y = xp(2);
z = xp(3);
vx = vp(1);
vy = vp(2);
vz = vp(3);

end