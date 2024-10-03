%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script starts with an orbit found via the
% 'commensurated_sigma_homotopy.m' script that exists at a RAAN value of
% zero.  Via a variable-stepsize natural parameter continuation process,
% this script increases this RAAN and reconverges for periodic solutions
% along the way.
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
RAAN = 0;
inc = deg2rad(5.145);

%% Definition of option structures
opts = odeset("RelTol", 1e-12, "AbsTol", 1e-12);

%% Load in the pre-existing sigma homotopy
% Load in the pre-existing data
load('..\saved data\generated\commensurate_homotopy_sigma.mat')

% Which dataset orbit is closest to our desired starting RAAN value
orbit = orbit(742);
orbit.RAAN = 0;
orbit.M_0 = 0;

% Solution structure to organize MS data on each propagated arc
int_results = rmfield(orbit.int_results, 'dxf_dsigmak');

% We know that 115 revs of our CR3BP SPO equal 240 * pi n.d. time units
number_arcs = 115; % This can however be anything (as long as M-S converges)

% How much elapsed time occurs on each arc?
delta_tau = orbit.TIP / number_arcs;

% Preallocate X (our design variable vector) with appropriate amt. of space
% We are doing full-state continuity only, with no augmentation
X = zeros([6 * number_arcs, 1]);

for k = 1:1:number_arcs
    X(6*k-5:6*k) = int_results(k).x_0k; % Fill 'X' with pre-existing data
end

%% Proceed through the multiple shooting method
% What value of sigma do we want an orbit at?
RAAN_target = 0;

% How many revs through our continuation process do we want to perform?
M_start = 1;
M_end = 100000;

% What NPC stepsize do we want to start at?  Changes adaptively.
ds = 1e-4;

% When is our N-R scheme successful?
constraint_tolerance = 1e-10;

% Do we want to plot every converged answer we find?
plot_converged_solns = false;

% Do we want to plot every trajectory propagated by the MS scheme?
plot_intermediates = false;

for M = M_start:M_end
    q = 1; % How many revs thru the N-R scheme have we performed?
    
    % Dummy initializations for the constraint vec. to start the N-R process
    F_norm_srz = rand(1, 3);
    F = 1;

    while norm(F) > constraint_tolerance
        % Propagate the orbit for the sigma value in the 'X' vector
        parfor k = 1:number_arcs
            % Shift the starting ArgLat. and relative angle B with arc's epoch
            starting_M = M0 + int_results(k).epoch;
            starting_B = RAAN_target - sqrt((muS + 1)/aE^3) * int_results(k).epoch;
    
            % Define initial conditions and anonymous func. for integration
            sv_k = [int_results(k).x_0k; reshape(eye(10), [100, 1])];
    
            % 'M0' and 'RAAN' for this sim start shifted wrt set values above
            ode_func = @(t, y) bcir4bp_stm(t, y, ...
                                        earth_moon_massparam = mu, ...
                                        sun_sgp_nondim = muS, ...
                                        sun_effect_slider = 1.0, ...
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
            % into it's constitutent basal STM
            augmented_STMf = reshape(ssk.y(7:106, end), [10, 10]);
    
            STMf = augmented_STMf(1:6, 1:6);
    
            % Pack up all of the results into the solution structure
            int_results(k).x_fk = x_fk;
            int_results(k).stm_fk = STMf;
        end
    
        % From the integrated information, construct 'F' and 'DF' appropriately
        F = zeros([6 * number_arcs, 1]);
        DF = zeros(6 * number_arcs);
    
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
        end
    
        % Diagnostic output of |F|
        fprintf("Norm of constraint vector:  %1.4e\n", norm(F))
    
        % Solve this linear system the 'best' way possible
        [delta_x, step_residual] = robust_delta_X_DG_G(DF, F);

        % Compute NPC stepsize adaptation factors
        if q == 1 % If this is the first rev through the N-R process
            first_SL = norm(delta_x); % How big was our step?
        end

        if q == 2 % If this is the second rev through the N-R process
            contr_rate = norm(delta_x) / first_SL; % Stepsize changing fast?
        end

        % What angle does current nullspace make with continuation nullspace?
        % NOTE: This stays zero because we're using NPC and not PAC
        step_angle = 0;
    
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
    
    if plot_converged_solns == true
        for k = 1:1:number_arcs
            % Shift the starting ArgLat. and relative angle B with arc's epoch
            starting_M = M0 + int_results(k).epoch;
            starting_B = RAAN_target - sqrt((muS + 1)/aE^3) * int_results(k).epoch;
    
            % Define initial conditions and anonymous func. for integration
            sv_k = int_results(k).x_0k;
    
            % 'M0' and 'RAAN' for this sim start shifted wrt set values above
            ode_func = @(t, y) bcir4bp_stm(t, y, ...
                                            earth_moon_massparam = mu, ...
                                            sun_sgp_nondim = muS, ...
                                            sun_effect_slider = 1.0, ...
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

    % Write converged orbit to organizing structure
    orbit(M).int_results = int_results;
    orbit(M).number_arcs = number_arcs;
    orbit(M).corrector = strcat('Fixed-time, full-rev patch point ', ...
                                'continuity only over a RAAN NPC');
    orbit(M).RAAN = RAAN_target;
    orbit(M).M_0 = 0;
    orbit(M).epoch = 0;
    orbit(M).sigma = 1.0;
    orbit(M).F_norm = norm(F);
    orbit(M).TIP = orbit(1).TIP;

    if M == 1
        orbit = orderfields(orbit);
        
        % Preallocate the structure for space, will be trimmed afterwards
        orbit(M_end).epoch = -1;
    end

    % Break out of the PAC process if we are finished
    if X(end) > 2 * pi
        break
    end

    %% Continue on with steplength adaptation and taking a PAC step
    % Nominal values for steplength adaptation, compared against actuals
    nominal_contr_rate = 1.3; % Larger => bigger ds
    nominal_first_SL = 2e-2;
    nominal_step_angle = 1e-1; % 0.1 radians is approximately 11 degrees

    % Add some hard-coded bounds to cap/lower bound the adaptation factor
    maximum_adaptation = 2;
    minimum_adaptation = 0.5;

    % Add some hard-coded bounds to keep the PAC scheme from moving fast/slow
    minimum_steplength = 1e-3;
    maximum_steplength = 0.1;

    % If the chosen member was converged on iter #1, we need to set a
    % dummy contraction rate for that iteration.  Choose one that allows
    % for the other convergence metrics to be applied instead.
    if ~exist("contr_rate", "var")
        contr_rate = first_SL / nominal_first_SL * nominal_contr_rate;
    end

    % When a bifurcation occurs, the step angle is closer to pi
    step_angle = min(step_angle, pi-step_angle);

    % Three factors can be used to adapt our steplength
    adapt_metrics = [sqrt(contr_rate/nominal_contr_rate), ...
                     sqrt(first_SL/nominal_first_SL), ...
                     step_angle / nominal_step_angle];

    largest_factor = max(adapt_metrics);

    % The most limiting of the three is chosen as our actual adaptation
    adapt_factor = max( min(largest_factor, maximum_adaptation), ...
                        minimum_adaptation);

    ds = min(max(ds/adapt_factor, minimum_steplength), maximum_steplength);

    % Print converged RAAN value and arclength for logging
    fprintf("RAAN = %5.4e | ds = %5.4e | M = %d | q = %d\n", RAAN_target, ds, M, q)

    % Update the RAAN target for the next step
    RAAN_target = RAAN_target + ds;
end