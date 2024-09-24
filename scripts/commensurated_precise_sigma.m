%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script starts with an orbit found via the
% 'commensurated_sigma_homotopy.m' script that exists near to a target
% value of sigma for which we'd like an orbit.  Using this starting
% orbit, the script corrects the initial conditions (and the value of
% sigma) with the intent of correcting the orbit and sigma to match the
% precise value of sigma specified
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

%% Load in the pre-existing sigma homotopy
% Load in the pre-existing data
load('..\saved data\generated\commensurate_homotopy_sigma.mat')

% Which dataset orbit is closest to our desired sigma value
M_close = 1;

% Solution structure to organize MS data on each propagated arc
int_results = struct('x_0k', {}, ...
                     'epoch', {}, ...
                     'int_time', {}, ...
                     'x_fk', {}, ...
                     'dxf_dsigmak', {}, ...
                     'stm_fk', {});

% We know that 115 revs of our SPO equal 240 * pi n.d. time units
number_arcs = 115; % This can however be anything (as long as M-S converges)

% How much elapsed time occurs on each arc?
delta_tau = T_repeat_nd / number_arcs;

% Preallocate X (our design variable vector) with appropriate amt. of space
X = zeros([6 * number_arcs + 1, 1]);

for k = 1:1:number_arcs
    % Overwrite the 'X' vector with pre-existing data
    X(6*k-5:6*k) = orbit(M_close).int_results(k).x_0k;

    % Seed the MS structure with proper values as well
    int_results(k).x_0k = orbit(M_close).int_results(k).x_0k;
    int_results(k).epoch = (k - 1) * delta_tau;
    int_results(k).int_time = delta_tau;
end

% Set the value of sigma properly too
X(end) = orbit(M_close).sigma;

%% Proceed through the multiple shooting method
% What value of sigma do we want an orbit at?
sigma_target = 0.0;

% When is our N-R scheme successful?
constraint_tolerance = 1e-12;

% Do we want to plot every converged answer we find?
plot_converged_solns = false;

% Do we want to plot every trajectory propagated by the MS scheme?
plot_intermediates = false;

q = 1; % How many revs thru the N-R scheme have we performed?

% Dummy initializations for the constraint vec. to start the N-R process
F_norm_srz = rand(1, 3);
F = 1;

while norm(F) > constraint_tolerance
    % Propagate the orbit for the sigma value in the 'X' vector
    parfor k = 1:number_arcs
        % Shift the starting ArgLat. and relative angle B with arc's epoch
        starting_M = M0 + int_results(k).epoch;
        starting_B = RAAN - sqrt((muS + 1)/aE^3) * int_results(k).epoch;

        % Define initial conditions and anonymous func. for integration
        sv_k = [int_results(k).x_0k; reshape(eye(7), [49, 1])];

        % 'M0' and 'RAAN' for this sim start shifted wrt set values above
        ode_func = @(t, y) bcir4bp_stm(t, y, ...
                                    earth_moon_massparam = mu, ...
                                    sun_sgp_nondim = muS, ...
                                    sun_effect_slider = X(end), ...
                                    moon_arglat_at_epoch = starting_M, ...
                                    moon_inclination = inc, ...
                                    moon_right_ascension = starting_B, ...
                                    earth_sma_nondim = aE, ...
                                    stm_enabled = true);

        % Return a solution structure corresponding to arc 'k'
        ssk = ode89(ode_func, [0, int_results(k).int_time], sv_k, opts);

        if plot_intermediates == true
            % Plot the unconverged MS arcs
            plot3(ssk.y(1, :), ssk.y(2, :), ssk.y(3, :), 'r')
        end

        % Unpack the final state of the propagation
        x_fk = ssk.y(1:6, end);

        % Unpack the augmented state transition matrix and decompose
        % into it's constitutent partial derivative matrices and vectors
        augmented_STMf = reshape(ssk.y(7:55, end), [7, 7]);

        STMf = augmented_STMf(1:6, 1:6);
        dxf_dsigma = augmented_STMf(1:6, 7);

        % Pack up all of the results into the solution structure
        int_results(k).x_fk = x_fk;
        int_results(k).stm_fk = STMf;
        int_results(k).dxf_dsigmak = dxf_dsigma;
    end

    % From the integrated information, construct 'F' and 'DF' appropriately
    F = zeros([6 * number_arcs + 1, 1]);
    DF = zeros(6 * number_arcs + 1);

    for k = 1:1:(number_arcs)
        index_range_a = (6 * k - 5):1:(6 * k);
        index_range_b = (6 * k + 1):1:(6 * k + 6);

        if k < number_arcs
            F(index_range_a) = int_results(k).x_fk - int_results(k+1).x_0k;

            % Tile the STM and an identity matrix in the right spots
            DF(index_range_a, index_range_a) = int_results(k).stm_fk;
            DF(index_range_a, index_range_b) = -eye(6);
        else
            F(index_range_a) = int_results(k).x_fk - int_results(1).x_0k;

            % The periodicity condition is enforced here, wraps to #1 IC
            DF(index_range_a, index_range_a) = int_results(k).stm_fk;
            DF(index_range_a, 1:6) = -eye(6);
        end

        % Add in the sigma partial derivative information
        DF(index_range_a, end) = int_results(k).dxf_dsigmak;
    end

    % Apply the constraint for a specific value of sigma
    F(end) = X(end) - sigma_target;
    DF(end, end) = 1;

    % Diagnostic output
    fprintf("Norm of constraint vector:  %1.4e\n", norm(F))

    % Solve this system the 'best' way possible
    [delta_x, step_residual] = robust_delta_X_DG_G(DF, F);

    % Prematurely break out of the N-R algorithm if constraint norm
    % is satisfied by our previous step (avoiding applying another dX)
    if norm(F) < constraint_tolerance
        break
    end

    % Modify our design variables
    X = X - delta_x;

    % Check that F is changing (e.g. dropping monotonically)
    % This is a moving average of the constraint vector's norm
    F_norm_srz = [F_norm_srz(2:3), norm(F)];
    average_change = mean(abs(diff(F_norm_srz))) / norm(F);

    if average_change < 1e-4
        error("F is not falling monotonically - breaking...")
    end

    % Deal the modified initial conditions back out to the struct
    for k = 1:1:number_arcs
        int_results(k).x_0k = X(6*k-5:6*k);
    end

    q = q + 1;

    if q > 50
        error("N-R scheme took more than 50 revs - breaking...")
    end
end

%% Continue with PAC procedure, first doing some preliminaries
% Plot the converged MS solution
if plot_converged_solns == true
    for k = 1:1:number_arcs
        % Shift the starting ArgLat. and relative angle B with arc's epoch
        starting_M = M0 + int_results(k).epoch;
        starting_B = RAAN - sqrt((muS + 1)/aE^3) * int_results(k).epoch;

        % Define initial conditions and anonymous func. for integration
        sv_k = int_results(k).x_0k;

        % 'M0' and 'RAAN' for this sim start shifted wrt set values above
        ode_func = @(t, y) bcir4bp_stm(t, y, ...
                                    earth_moon_massparam = mu, ...
                                    sun_sgp_nondim = muS, ...
                                    sun_effect_slider = X(end), ...
                                    moon_arglat_at_epoch = starting_M, ...
                                    moon_inclination = inc, ...
                                    moon_right_ascension = starting_B, ...
                                    earth_sma_nondim = aE, ...
                                    stm_enabled = false);

        % Return a solution structure corresponding to arc 'k'
        ssk = ode45(ode_func, [0, int_results(k).int_time], sv_k, opts);

        % Plot the converged MS solution
        plot3(ssk.y(1, :), ssk.y(2, :), ssk.y(3, :), 'b')
    end
end

orbit(M_close).F_norm = norm(F);
orbit(M_close).sigma = sigma_target;
orbit(M_close).int_results = int_results;