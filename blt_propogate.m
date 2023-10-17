% Code by Andrew Binder

%% MATLAB Initialization and MATLAB Constants Definition
clear; clc; close('all')

addpath(genpath('utilities'));  % Add folder and subfolders
addpath(genpath('saved data'));  % Add folder and subfolders

% Sphere mesh points for plotting of the primaries
[X_sphere, Y_sphere, Z_sphere] = sphere(10);

%% Constants of the Problem
mu_moon = 4902.8005821478;  % All in km^3 / s^2
mu_earth = 398600.4415;
mu_sun = 132712440018;

l_star = 384400;  % Also equal to the SMA of the Moon
t_star = sqrt(l_star^3 / (mu_moon+mu_earth)); % Divide by 86400 to get time in days
v_star = l_star / t_star;

a_s = 149.598 * 10^6 / l_star;  % E-M Barycenter's SMA about the Sun, [nd]

mu = mu_moon / (mu_earth + mu_moon);  % mu value for the Earth-Moon system, [nd]
m_s = mu_sun / (mu_earth + mu_moon);  % Nondimensionalized solar mass

mu_oom = round(log10(mu)); mu_SF_known = 8;  % We know 8 Earth sigfigs, 9 Moon sigfigs, and 9 Sun sigfigs

% mu rounded to appropriate significance
mu = round(mu*10^(mu_SF_known-mu_oom-1)) / 10^(mu_SF_known-mu_oom-1);

%% Definition of option structures

% Integrate to a tight tolerance - long-term sims
opts = odeset("RelTol", 3e-14, "AbsTol", 3e-14, "MaxStep", 0.03);

%% Initialize a trajectory figure
% Orbital figure initialize - lib. points and E/M
orbitviews.shalos = figure(); shalos_plots = [];
axis equal; grid on; hold on;

shalos_plots.earth = surf(X_sphere * 6378/l_star - mu, Y_sphere * 6378/l_star, Z_sphere * 6378/l_star, ...
            'EdgeColor', '#0047AB', 'FaceColor', '#6495ED');
shalos_plots.moon = surf(X_sphere * 1738/l_star + 1 - mu, Y_sphere * 1738/l_star, Z_sphere * 1738/l_star, ...
            'EdgeColor', '#696969', 'FaceColor', '#848884');

xlabel("x-distance from barycenter [n.d.]")
ylabel("y-distance from barycenter [n.d.]")
zlabel("z-distance from barycenter [n.d.]")
% End orbital figure initialize

%% BLT Propogation

% Load in one of the nice ballistic lunar transfers
load("bcr4bp\saved data\generated\periapsis_maps\selected_BLT_states\R21S2_BLT_1.mat")

% Properly choose 'epoch' Sun angle (RPO insertion @ t = 0)
th_S0 = blt_prop_time * (sqrt((m_s + 1)/a_s^3) - 1);

% Integrate from the Earth until RPO insertion
sol_struct = ode89(@(t,y) x_bar_dot(t, y, mu, m_s, a_s, th_S0), [0, -blt_prop_time], blt_IC, opts);

% Plot this full simulation to visualize 'settling' period
settle = plot3(sol_struct.y(1, :), sol_struct.y(2, :), sol_struct.y(3, :), 'b')

% T = 15.8 for BLT1, 8.5 for BLT2
% Also produce a truncated sim from the Earth until lunar flyby
sol_struct = ode89(@(t,y) x_bar_dot(t, y, mu, m_s, a_s, th_S0), [0, 8.5], blt_IC, opts);

% Plot traj. until lunar flyby + lunar close approach pt.
blt = plot3(sol_struct.y(1, :), sol_struct.y(2, :), sol_struct.y(3, :), 'g')
flyby = scatter3(sol_struct.y(1, end), sol_struct.y(2, end), sol_struct.y(3, end), 'filled', 'g')

% Also mark the RPO injection point
inject = scatter3(blt_FC(1), blt_FC(2), blt_FC(3), 'filled', 'ks')

% Simulate the 'RPO' (unconverged in BCR4BP) for vis. purposes
blt_FC(4:6) = (1 - 0.1 / norm(blt_FC(4:6) / v_star)) * blt_FC(4:6)
sol_struct = ode89(@(t,y) x_bar_dot(t, y, mu, m_s, a_s, 0), [0, 4 * pi], blt_FC, opts);

% Plot the 'RPO'
orbit = plot3(sol_struct.y(1, :), sol_struct.y(2, :), sol_struct.y(3, :), 'r')