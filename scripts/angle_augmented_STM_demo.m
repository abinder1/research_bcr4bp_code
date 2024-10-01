%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is a demonstration of the Lunar orientation angle state
% transition matrix derivatives.  Each of the three new derivatives
% (dxf/dM0, dxf/d\Omega, and dxf/di) are calculated and checked against
% finite difference approximations of each derivative
%
% Author:  Andrew Binder (2024)
%
% Inputs: None
% Ouptuts: None
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

aE = aE_dim / l_star; % Nondimensionalize the Earth-Sun distance

% Lunar orientation relative to ecliptic inertial frame at epoch
M0 = 0;
inc = deg2rad(5.145);
RAAN = 0;

%% Definition of option structures
opts = odeset("RelTol", 1e-12, "AbsTol", 1e-12);

%% Propagate a sample reference trajectory in the BCIR4BP

% Load in L4 SPOs, our sample is entry #48
load('..\saved data\generated\l4_short_period.mat')

%% Pull a sample and integrate for the new derivatives
% Choose the commensurate orbit sample from the dataset (#48 is our
% newly-constructed commensurate SPO member)
orbit = l4_short_period(48);

% Propagate the orbit at full sun-strength to get the new partials
ode_func = @(t, y) bcir4bp_angles_stm(t, y, ...
                                      earth_moon_massparam = mu, ...
                                      sun_sgp_nondim = muS, ...
                                      sun_effect_slider = 1.0, ...
                                      moon_arglat_at_epoch = M0, ...
                                      moon_inclination = inc, ...
                                      moon_right_ascension = RAAN, ...
                                      earth_sma_nondim = aE, ...
                                      stm_enabled = true);

T_prop = 5; % Reasonably long propagation time
sv_0 = [orbit.ic; reshape(eye(10), [100, 1])];
sol_struct = ode45(ode_func, [0, T_prop], sv_0, opts);

% Unpack the final state of the propagation
x_fk = sol_struct.y(1:6, end);

% Unpack the augmented state transition matrix and decompose
% into it's constitutent partial derivative matrices and vectors
augmented_STMf = reshape(sol_struct.y(7:106, end), [10, 10]);

STMf = augmented_STMf(1:6, 1:6);
dxf_dsigma = augmented_STMf(1:6, 7);
dxf_dM0 = augmented_STMf(1:6, 8);
dxf_dOmega = augmented_STMf(1:6, 9);
dxf_dinc = augmented_STMf(1:6, 10);

%% Check the partials with finite differencing
delta = 1e-6; % The size of the finite difference step

% First check the M0 derivatives
% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_angles_stm(t, y, ...
                                      earth_moon_massparam = mu, ...
                                      sun_sgp_nondim = muS, ...
                                      sun_effect_slider = 1.0, ...
                                      moon_arglat_at_epoch = M0+delta, ...
                                      moon_inclination = inc, ...
                                      moon_right_ascension = RAAN, ...
                                      earth_sma_nondim = aE, ...
                                      stm_enabled = false);

sol_struct_perturbed = ode45(ode_func, [0, T_prop], orbit.ic, opts);

% Unpack the final state of the propagation
x_fk_M0_perturbed = sol_struct_perturbed.y(1:6, end);

dxf_dM0_fd = (x_fk_M0_perturbed - x_fk) / delta
dxf_dM0

% Next check the \Omega derivatives
% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_angles_stm(t, y, ...
                                      earth_moon_massparam = mu, ...
                                      sun_sgp_nondim = muS, ...
                                      sun_effect_slider = 1.0, ...
                                      moon_arglat_at_epoch = M0, ...
                                      moon_inclination = inc, ...
                                      moon_right_ascension = RAAN+delta, ...
                                      earth_sma_nondim = aE, ...
                                      stm_enabled = false);

sol_struct_perturbed = ode45(ode_func, [0, T_prop], orbit.ic, opts);

% Unpack the final state of the propagation
x_fk_M0_perturbed = sol_struct_perturbed.y(1:6, end);

dxf_dOm_fd = (x_fk_M0_perturbed - x_fk) / delta
dxf_dOmega

% And finally check the inclination derivatives
% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_angles_stm(t, y, ...
                                      earth_moon_massparam = mu, ...
                                      sun_sgp_nondim = muS, ...
                                      sun_effect_slider = 1.0, ...
                                      moon_arglat_at_epoch = M0, ...
                                      moon_inclination = inc+delta, ...
                                      moon_right_ascension = RAAN, ...
                                      earth_sma_nondim = aE, ...
                                      stm_enabled = false);

sol_struct_perturbed = ode45(ode_func, [0, T_prop], orbit.ic, opts);

% Unpack the final state of the propagation
x_fk_M0_perturbed = sol_struct_perturbed.y(1:6, end);

dxf_dIn_fd = (x_fk_M0_perturbed - x_fk) / delta
dxf_dinc