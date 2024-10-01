%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script finds an orbit with the commensurate period correct for
% use in the commensurate-d model of the BCIR4BP.  Using a value of
% sigma = 0.0, the BCIR4BP and it's associated STM are effectively
% simulating CR3BP dynamics and the CR3BP STM.  The shooting process is
% done here to demonstrate the equivalence of the BCIR4BP model to the
% CR3BP model when sigma = 0.0.
%
% Author:  Andrew Binder (2024)
%
% Inputs: None
% Outputs: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB Initialization and MATLAB Constants Definition
clear;
clc;
close('all');

addpath(genpath('..\utilities')); % Add folder and subfolders
addpath(genpath('..\saved data')); % Add folder and subfolders

%% Constants of the Problem
mu_moon = 4902.8005821478; % All in km^3 / s^2
mu_earth = 398600.4415;
mu_sun = 132712440018;

aE_dim = 149598023; % Earth's SMA about the Sun

mu = mu_moon / (mu_earth + mu_moon); % mu value for the Earth-Moon system
mu_oom = round(log10(mu));
mu_SF_known = 8; % We know 8 Earth sigfigs, 9 Moon sigfigs

% mu rounded to appropriate significance
mu = round(mu*10^(mu_SF_known - mu_oom - 1)) / 10^(mu_SF_known - mu_oom - 1);

muS = mu_sun / (mu_earth + mu_moon); % Sun's n.d. SGP [L^3 / T^2]

% Commensurable char. length has a t* that is best rational approx. to a tol.
N = 3;
D = 40;

% Approximately equal to the average SMA of the Moon about EM Bary
l_star = ((N / D)^2 / (muS + 1))^(1 / 3) * aE_dim;

% l_star = 384400;  % Equal to the average SMA of the Moon about EM Bary
t_star = sqrt(l_star^3/(mu_moon + mu_earth)); % Div. by 86400 > time in days
v_star = l_star / t_star; % Characteristic speed, in km/s

aE = aE_dim / l_star;  % Nondimensionalize the Earth-Sun distance

% Lunar orientation relative to ecliptic inertial frame at epoch
M0 = 0;
inc = deg2rad(5.145);
RAAN = 0;

%% Definition of option structures
opts = odeset("RelTol", 1e-12, "AbsTol", 1e-12);

%% Propagate a sample reference trajectory in the BCIR4BP

% Load in L4 SPOs for demonstration of commensurated dynamics
load('..\saved data\generated\l4_short_period.mat')

% Find the Q value such that 2 * pi * D / Q is as close as possible to
% the period of 'SPO1' ~= index #44 of my dataset
Q = 115;

% An SPO with this synodic period repeats about every nine sidereal years
TIP_commens = 2 * (3* D) * pi / Q; % Considerably reduced from 'commens_demo.m'

% M(\tau), B(\tau), and the chosen PO all repeat at this nondimensional period
T_repeat_nd = 240 * pi; % Considerably reduced from 'commens_demo.m'

% Choose the commensurate orbit's proxy sample from the dataset
orbit = l4_short_period(47);

%% Try to solve for the IC with the right commensurate period
% Set the Sun's effect to no-strength for the single-shooting process
ode_func = @(t, y) bcir4bp_stm(t, y, ...
                                      earth_moon_massparam = mu, ...
                                      sun_sgp_nondim = muS, ...
                                      sun_effect_slider = 0.0, ...
                                      moon_arglat_at_epoch = M0, ...
                                      moon_inclination = inc, ...
                                      moon_right_ascension = RAAN, ...
                                      earth_sma_nondim = aE, ...
                                      stm_enabled = true);

% Initialize the guess with a good nearby IC
X = orbit.ic;

% Dummy initialization of the F vector to get the loop going
F = 1;
k = 1;

figure(1); axis equal; hold on; grid on;

% Find the CR3BP IC that is periodic at the right period
while norm(F) > 1e-12
    sv0 = [X; reshape(eye(10), [100, 1])];
    bcir_ss = ode45(ode_func, [0, TIP_commens], sv0, opts);

    plot3(bcir_ss.y(1, :), bcir_ss.y(2, :), bcir_ss.y(3, :))
    
    statef = bcir_ss.y(1:6, end);
    augmented_STMf = reshape(bcir_ss.y(7:106, end), [10, 10]);

    STMf = augmented_STMf(1:6, 1:6);
    dxf_dsigma = augmented_STMf(1:6, 7);
    
    F = statef - X;
    DF = STMf - eye(6);
    
    deltaX = robust_delta_X_DG_G(DF, F);
    
    X = X - deltaX;

    k = k + 1;
end

% Construct new entry for insertion into the L4 SPO data structure
new_l4_spo_entry.TIP = TIP_commens;
new_l4_spo_entry.ic = X;
new_l4_spo_entry.k = k;
new_l4_spo_entry.method = '09/22/24 Binder L4 Sh. Pd. Full-State SS Targeter';
new_l4_spo_entry.mono = STMf;

% New entry is at element #48 in the data structure
l4_short_period = [l4_short_period(1:47), ...
                   new_l4_spo_entry, ...
                   l4_short_period(49:end)];

% Visualize the CR3BP's propogation of our commensurate IC and mark its initial
% condition and its final condition
figure(1);

hold on;
axis equal;
grid on;

plot3(bcir_ss.y(1, :), bcir_ss.y(2, :), bcir_ss.y(3, :), 'black')
scatter3(bcir_ss.y(1, end), bcir_ss.y(2, end), bcir_ss.y(3, end), 'rs')
scatter3(bcir_ss.y(1, 1), bcir_ss.y(2, 1), bcir_ss.y(3, 1), 'ks')

xlabel("x-distance from EM barycenter [nd]")
ylabel("y-distance from EM barycenter [nd]")
zlabel("z-distance from EM barycenter [nd]")

title("Shooting process to find commensurate SPO")