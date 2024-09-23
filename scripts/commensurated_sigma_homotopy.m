%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script starts with an orbit from the CR3BP, the commensurate L4
% SPO, and rev stacks it at a value of sigma = 0.0.  Using this
% rev-stacked orbit and using the partial derivatives with respect to
% sigma, this script follows a pseudoarclength continuation process to
% continue this orbit until a value of sigma = 1.0.  This converged
% orbit is then saved to file, and plotted for visualization's sake.
% This process is followed at the correct lunar inclination of 5.145
% degrees, but at values of M0 and RAAN of zero - continuation over
% these variables will be done later (once a sigma = 1.0 orbit has been
% found).
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

% M(\tau), B(\tau), and the chosen PO all repeat at this nondimensional period
T_repeat_nd = 240 * pi; % Considerably reduced from 'commens_demo.m'

%% Build out the right initial X vector via rev stacking
% Choose the commensurate orbit sample from the dataset
orbit = l4_short_period(48);

% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_stm(t, y, ...
                               earth_moon_massparam = mu, ...
                               sun_sgp_nondim = muS, ...
                               sun_effect_slider = 0.0, ...
                               moon_arglat_at_epoch = M0, ...
                               moon_inclination = inc, ...
                               moon_right_ascension = RAAN, ...
                               earth_sma_nondim = aE, ...
                               stm_enabled = false);

single_rev_sol_struct = ode45(ode_func, [0, orbit.TIP], orbit.ic, opts);

figure; hold on; axis equal; grid on;
plot3(single_rev_sol_struct.y(1,:), ...
      single_rev_sol_struct.y(2,:), ...
      single_rev_sol_struct.y(3,:), 'k')

% Solution structure to organize MS data on each propagated arc
int_results = struct('x_0k', {}, ...
                     'epoch', {}, ...
                     'int_time', {}, ...
                     'x_fk', {}, ...
                     'dxf_dsigmak', {}, ...
                     'stm_fk', {});

% We know that 115 revs of our SPO equal 240 * pi n.d. time units
number_arcs = 115;

% How much elapsed time occurs on each arc?
delta_tau = T_repeat_nd / number_arcs;

% Preallocate X (our design variable vector) with appropriate amt of space
X = zeros([6 * number_arcs + 1, 1]);

% Sample the singly-propagated trajectory for many 'revs'
for k = 1:1:number_arcs
    epoch_time = (k-1) * delta_tau;
    time_on_orbit = mod(epoch_time, orbit.TIP);

    % Populate X six indices at a time with samples from trajectory
    starting_state = deval(single_rev_sol_struct, time_on_orbit, 1:6);

    % Save the epoch, initial state, and integration time on each arc
    int_results(k).epoch = epoch_time;
    int_results(k).x_0k = starting_state;
    int_results(k).int_time = delta_tau;

    X(6*k-5:6*k) = starting_state;
end

% Choose an initial sigma value of zero
X(end) = 0;

%% Construct the multiple-shooting algorithm for the homotopy process