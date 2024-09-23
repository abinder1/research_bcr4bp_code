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
G_norm_srz = rand(1,3); G = 1;  % Dummy initializations for the constraint vec.

% How many revs through our continuation process do we want to perform?
M_start = 1;
M_end = 1;

% Which col. of the singular value decomposition do we consider our nullspace?
tan_vec_select = length(X);

% What pseudoarclength stepsize do we want to start at?  Changes adaptively
ds = 1e-4;

% When is our N-R scheme successful?
constraint_tolerance = 1e-12;

% Do we want to plot every converged answer we find?
plot_converged_solns = true;

% Do we want to plot every trajectory propagated by the MS scheme?
plot_intermediates = true;

for M = M_start:M_end
    q = 1;  % How many revs thru the N-R scheme have we performed?

    while norm(G) > constraint_tolerance
        % Propagate the orbit for the sigma value in the 'X' vector
        parfor k = 1:number_arcs
            % Shift the starting Arg. Lat. and relative angle B with arc's epoch
            starting_M = M0 + int_results(k).epoch;
            starting_B = RAAN - sqrt((muS + 1)/aE^3) * int_results(k).epoch;
    
            % Define initial conditions and anonymous func. for integration
            sv_k = [int_results(k).x_0k; reshape(eye(7), [49,1])];
    
            % 'M0' and 'RAAN' for this sim start shifted wrt their normal values
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
            ssk = ode45(ode_func, [0, int_results(k).int_time], sv_k, opts);

            if plot_intermediates == true
                % Plot the unconverged MS arcs
                plot3(ssk.y(1,:), ssk.y(2,:), ssk.y(3, :), 'r')
            end
    
            % Unpack the final state of the propagation
            x_fk = ssk.y(1:6, end);
    
            % Unpack the augmented state transition matrix and decompose
            % into it's constitutent partial derivative matrices and vectors
            augmented_STMf = reshape(ssk.y(7:55, end), [7 7]);
    
            STMf = augmented_STMf(1:6, 1:6);
            dxf_dsigma = augmented_STMf(1:6, 7);
    
            % Pack up all of the results into the solution structure
            int_results(k).x_fk = x_fk;
            int_results(k).stm_fk = STMf;
            int_results(k).dxf_dsigmak = dxf_dsigma;
        end
    
        % From the integrated information, construct 'F' and 'DF' appropriately
        F = zeros([6 * number_arcs, 1]);
        DF = zeros([6 * number_arcs, 6 * number_arcs + 1]);
    
        for k = 1:1:(number_arcs)
            index_range_a = (6*k-5):1:(6*k);
            index_range_b = (6*k+1):1:(6*k+6);
    
            if k < number_arcs
                F(index_range_a) = int_results(k).x_fk - int_results(k + 1).x_0k;
    
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
    
        % Find the singular values and nullspace vector for pseudoarclength
        [~, S, V] = svd(DF);
    
        if M > M_start
            % Determine if my current 'X' vector lies along the nullspace
            % at a distance of 'ds' from my last converged answer 'last_X'
            F_pseudoarc = dot(X - last_X, step_tangent_vector) - ds;
        else
            % If we are on the first PAC rev and first iter through it,
            % we want to set our 'tangent_vector' or the unit vector
            % pointing into our nullspace
            if q == 1  
                tangent_vector = V(:, tan_vec_select)';

                % 'step_tangent_vector' is the unit vector oriented in
                % the correct pointing direction (pos or neg)
                step_tangent_vector = tangent_vector;
            end
    
            % There is no previous converged answer to compare against,
            % so for M = M_start this constraint is zero
            F_pseudoarc = 0;
        end

        % Make the DF matrix square - the nullspace is linearly indep.
        DG = [DF; step_tangent_vector];
        G = [F; F_pseudoarc];  % Append the step distance constraint

        % Diagnostic output
        fprintf("Norm of constraint vector:  %1.4e\n", norm(G))

        % Prematurely break out of the N-R algorithm if constraint norm
        % is satisfied by our previous step (avoiding applying another dX)
        if norm(G) < constraint_tolerance
            break
        end

        % Solve this system the 'best' way possible
        [delta_x, step_residual] = robust_delta_X_DG_G(DG, G);

        % Compute pseudo-arclength stepsize adaptation factors
        if q == 1  % If this is the first rev through the N-R process
            first_SL = norm(delta_x);  % How big was our step?
        end
    
        if q == 2  % If this is the second rev through the N-R process
            contr_rate = norm(delta_x) / first_SL;  % Stepsize changing fast?
        end
    
        if M > M_start  % Is this the first rev through PAC process?  If not:
            dot_prod = round(dot(step_tangent_vector, signum * tangent_vector), 13);
        else
            dot_prod = round(dot(tangent_vector, tangent_vector), 13);
        end
        
        % What angle does current nullspace make with continuation nullspace?
        step_angle = acos(dot_prod);

        % Modify our design variables
        X = X - delta_x;

        % Check that G is changing (e.g. dropping monotonically)
        % This is a moving average of the constraint vector's norm
        G_norm_srz(1) = [];  G_norm_srz = [G_norm_srz norm(G)];
        average_change = mean(abs(diff(G_norm_srz))) / norm(G);
    
        if average_change < 1e-4
            error("G is not falling monotonically - breaking...")
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

    %% Continue with PAC procedure once a soln is found
    % Plot the converged MS solution
    if plot_converged_solns == true
        for k = 1:1:number_arcs
            % Shift the starting Arg. Lat. and relative angle B with arc's epoch
            starting_M = M0 + int_results(k).epoch;
            starting_B = RAAN - sqrt((muS + 1)/aE^3) * int_results(k).epoch;
    
            % Define initial conditions and anonymous func. for integration
            sv_k = int_results(k).x_0k;
    
            % 'M0' and 'RAAN' for this sim start shifted wrt their normal values
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
            plot3(ssk.y(1,:), ssk.y(2,:), ssk.y(3, :), 'b')
        end
    end
end