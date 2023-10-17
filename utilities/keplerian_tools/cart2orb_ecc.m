function [a, ecc, inc, w, R, f] = cart2orb_ecc(x, y, z, vx, vy, vz)

%%  Coordinate transformation from cartesian coordinates to orbital elements.
%
%   Synopsis: this program transforms cartesian elements in the ECI frame
%             of reference into orbital elements (a, e, inc, w, R, f).
%
%   INPUT
%       x, y, z: cartesian position of the satellite (km)
%       vx, vy, vz: cartesian velocity of the satellite (km/s)
%
%   OUTPUT
%       a: semimajor axis (km)
%       ecc: eccentricity
%       inc: inclination (deg)
%       w: argument of perigee (deg)
%       R: right ascension of the ascending node (deg)
%       f: mean anomaly (deg)
%
%   Authors and Updates
%       David Arnas 
%       26 January 2022

%% Earth gravitational constant
mu = 398600.4418;

%% Radial distance
r = sqrt(x^2 + y^2 + z^2);

%% Angular moment
h = zeros(1,3);
h(1) = y*vz - z*vy;
h(2) = vx*z - x*vz;
h(3) = x*vy - y*vx;
mh = sqrt(h(1)^2 + h(2)^2 + h(3)^2);
uh = h/mh;

%% Eccentricity
ev = zeros(1,3);
ev(1) = (vy*h(3) - vz*h(2))/mu - x/r;
ev(2) = (vz*h(1) - vx*h(3))/mu - y/r;
ev(3) = (vx*h(2) - vy*h(1))/mu - z/r;
ecc = sqrt(ev(1)^2 + ev(2)^2 + ev(3)^2);
ne = ev/ecc;

%% Semi-major axis
a = (mh^2/mu)/(1 - ecc^2);

%% Inclination
inc = acos(h(3)/mh);

%% Line of nodes
n = zeros(1,3);
n(1) = -uh(2);
n(2) = uh(1);
n(3) = 0;
mn = sqrt(n(1)^2 + n(2)^2 + n(3)^2);
n = n /mn;

%% Right ascension of the ascending node
R = atan2(uh(1),-uh(2));

%% Argument of perigee
cosw = (ne(1)*n(1) + ne(2)*n(2));
sinw = ((n(2)*ne(3)-n(3)*ne(2))*uh(1) + (n(3)*ne(1)-n(1)*ne(3))*uh(2) + ...
    (n(1)*ne(2)-n(2)*ne(1))*uh(3));
w = atan2(sinw,cosw);

%% True anomaly
cosf = (ne(1)*x + ne(2)*y + ne(3)*z)/r;
sinf = ((ne(2)*z-ne(3)*y)*uh(1) + (ne(3)*x-ne(1)*z)*uh(2) + ...
    (ne(1)*y-ne(2)*x)*uh(3))/r;
f = atan2(sinf,cosf);

%% Tranformation to degrees
inc = inc*180/pi;
R = R*180/pi;
w = w*180/pi;
f = f*180/pi;

end