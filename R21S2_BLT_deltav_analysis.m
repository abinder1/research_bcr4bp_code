% Code by Andrew Binder

%% MATLAB Initialization and MATLAB Constants Definition
clear; clc; close('all');

addpath(genpath('utilities'));  % Add folder and subfolders
addpath(genpath('saved data'));  % Add folder and subfolders

%% Permit usage of NASA NAIF's MICE Toolkit
addpath(fullfile(pwd, './mice/src/mice/'))
addpath(fullfile(pwd, './mice/lib/'))

cspice_kclear;  % Clear loaded ephemeris kernels

%% Ephemeris file loading and access

% Get good absolute paths to the datafiles
planet_datafile = fullfile(pwd, './saved data/NASA NAIF/DE405AllPlanets.bsp');
leapseconds_datafile = fullfile(pwd, './saved data/NASA NAIF/naif0011.tls');
rotframe_datafile = fullfile(pwd, './saved data/NASA NAIF/frame kernels/EMBARY_ROTATING.tk');

% Load in the appropriate datafiles to MICE
cspice_furnsh( planet_datafile )  % Gives Earth, Moon positions/velocities
cspice_furnsh( leapseconds_datafile )  % Gives millisecond timing data
cspice_furnsh( rotframe_datafile )  % Gives millisecond timing data

%% Constants of the Problem
mu_moon = 4902.8005821478;  % Both in km^3 / s^2
mu_earth = 398600.4415;

l_star = 384400;  % Also equal to the average SMA of the Moon
t_star = sqrt(l_star^3 / (mu_moon+mu_earth)); % Divide by 86400 to get time in days
v_star = l_star / t_star;  % In km/s

mu = mu_moon / (mu_earth + mu_moon);  % mu value for the Earth-Moon system
mu_oom = round(log10(mu)); mu_SF_known = 8;  % We know 8 Earth sigfigs, 9 Moon sigfigs

% mu rounded to appropriate significance
mu = round(mu*10^(mu_SF_known-mu_oom-1)) / 10^(mu_SF_known-mu_oom-1);

%% Get POI location in J2000

% Load in the selected ballistic lunar transfer IC file
load('saved data\generated\periapsis_maps\selected_BLT_states\R21S2_BLT_2.mat')

% Set 't = 0' epoch to be this date - also the date of the next solar eclipse
epoch_burn = cspice_str2et('April 08, 2024 03:00 PM EDT');

% Figure out the date we leave LEO - ~100 days prior to epoch
epoch_depart = epoch_burn + blt_prop_time * 86400;

% Dimensionalized synodic state at Earth close approach
dim_depart_state = [blt_IC(1:3) * l_star; blt_IC(4:6) * v_star];

% This returns a pure rotation matrix, like the Q matrix from my notes on coordinate transforms
Q_mat = cspice_sxform('EM_BARYCENTRIC_ROT', 'J2000', epoch_depart );  % Translates from dim. synodic to J2000

% Now we can reference our new frame and get coordinate transformations
[syn_earth_state, ~] = cspice_spkezr( 'Earth', ...
                                      epoch_depart, ...
                                      'EM_BARYCENTRIC_ROT', ...
                                      'NONE', ...
                                      'EARTH-MOON BARYCENTER' );

% Rebase the synodic state and transform: from the barycenter to Earth-centered
MBJ2000_state = Q_mat * (dim_depart_state - syn_earth_state);

% What is the J2000 radial distance and altitude?
rad = norm(MBJ2000_state(1:3));  alt = rad - 6378;

% What is the satellite's speed?
V2 = norm(MBJ2000_state(4:6));

% Calculate the delta-v required for injection from a circular orbit 
V1 = sqrt(mu_earth / rad);  dV_POI = V2 - V1;

% Calculating Earth departure C3 is a helpful comparison metric with CAPSTONE
C3 = V2^2 - 2 * mu_earth / rad;

% Convert the hypothetical circular orbit J2000 state into Keplerian elems
[SMA, ecc, inc, AOP, RAAN, f] = cart2orb_ecc(   MBJ2000_state(1), ...
                                                MBJ2000_state(2), ...
                                                MBJ2000_state(3), ...
                                                MBJ2000_state(4) * V1 / V2, ...
                                                MBJ2000_state(5) * V1 / V2, ...
                                                MBJ2000_state(6) * V1 / V2      );

kep_c = [SMA; ecc; inc; AOP; RAAN; f]  % Compose into a vector for easy ref

% Both BLT's can be launched from Vandenberg AFB
% BLT #1 alt:  337.9 km     |     BLT #2 alt:  98.4 km
% BLT #1 inc:  111.2 deg    |     BLT #2 inc:  133.3 deg
% BLT #2 Vc:  7.7040 km/s   |     BLT #2 Vc:  7.8452 km/s
% BLT #1 dV:  3.1771 + 0.1  |     BLT #2 dV:  3.1763 + 0.1

%% "J2000" Definition:
% This ECI frame is centered at the Earth and has the Earth's velocity.
% According to the source below, most data products ref.d to J2000 are ICRF
% 
% X-axis:  The intersetion of the equatorial and ecliptic planes
% Y-axis:  Completes the triad
% Z-axis:  Normal to the mean equator of date at epoch J2000 TBD
%
%           - "J2000" in MICE is equivalent to GMAT's "EarthMJ2000Eq"
%
% Source:  https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/17_frames_and_coordinate_systems.pdf