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

N_m = 67;
N_e = 5;

% M(\tau), B(\tau), and the chosen PO all repeat at this nondimensional period
T_r = 2 * pi * N_m;

% Size of the STM (N x N)
N = 7;

% Define the initial lunar/solar angles for the homotopy
B0 = 1.5;
TH0 = 4.0;
sigma = 0;

%% Definition of option structures
opts = odeset("RelTol", 1e-13, "AbsTol", 1e-13);

%% Build out the right initial X vector via rev stacking
% Load in L4 SPOs, our sample is entry #48
load('..\saved data\generated\l4_short_period.mat')

% Choose the commensurate orbit sample from the dataset (#48 is our
% newly-constructed SPO member)
orbit = l4_short_period(131);

% Create a configuration for the simulation
sim_config = create_sim_config( sigma, ...  % \sigma
                                B0, ...  % B_0
                                TH0, ...  % \theta_0
                                N_m, ...  % N_m
                                N_e, ...  % N_e
                                0.0  );   % \tau_e

% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_stm(t, y, sim_config);

sv_0 = [orbit.ic; reshape(eye(sim_config.N), [sim_config.N^2, 1])];
single_rev_sol_struct = ode45(ode_func, [0, orbit.TIP], sv_0, opts);

%% Initialize the multiple shooting arcs and the X vector
N_arcs = 172;
prop_time = T_r / N_arcs;

X = zeros([N_arcs * 6, 1]);

for k = 1:1:N_arcs
    sim_config_temp = sim_config;
    sim_config_temp.sigma = 0;

    % All arcs are propagated for an even amount of time
    sim_config_temp.tau_e = (k-1) * prop_time;
    sim_config_temp.prop_time = prop_time;

    % Figure out where in the orbit this absolute tau is and pull state
    on_orbit_time = mod(sim_config_temp.tau_e, orbit.TIP);
    arc_initial_state = deval(single_rev_sol_struct, on_orbit_time, 1:6);

    % Put the pulled state into the design variable vector
    X(6*k-5:6*k) = arc_initial_state;

    if k == 1  % Preallocation
        sim_config_arcs(N_arcs) = sim_config_temp;
    end

    % Save configuration for this arc into a struct list
    sim_config_arcs(k) = sim_config_temp;
end

clear sim_config_temp

%% Begin correction process

% ----- Script output/visualization configuration ----- %
plot_converged_orbits = 0;  % Plot every X orbits

% ----- Newton-Raphson configuration/initialization ----- %
q_max = 20;  % How many N-R steps is too many?
convergence_tolerance = 1e-6;  % For the majority of N-R steps, what do we consider converged?

% Give the 'F' vector a dummy initialization
F = zeros([6*N_arcs, 1]);  
F = F + 1;

% Pre-allocate this structure using a dummy entry at the end
% (overwritten later)
propagations(N_arcs).x0_k = [0; 0; 0; 0; 0; 0];

% ----- Adaptive sigma modification configuration ----- %
min_ds = 1e-7;  % What is the smallest change in sigma we'll tolerate
fail_counter = 0;  % Adaptive stepsize tries and fails sigma values
ds = 1e-2;  % What change in sigma do we want to initialize with?

% When generating new sigma values after failure, we generate values in the range [10^min, 10^max]
ds_exponent_max = -1;
ds_exponent_min = log10(min_ds);

% These are index bounds on the number of sigma changes we plan to do
% This is done to avoid an infinite loop
M_start = 1;  M_end = max(10000, round(1 / min_ds));

%% Run the adaptive stepsize convergence process, starting with N-R
for M = M_start:1:M_end
    for q = 1:1:q_max
        % Preallocate the DF matrix as zeros
        DF = zeros([6*N_arcs, 6*N_arcs]);

        % Modify all arc's sigma values, and set initial conditions from X
        for k = 1:1:N_arcs
            sim_config_arcs(k).sigma = sigma;
            propagations(k).x0_k = X(6*k-5:6*k);
        end
    
        % Run integrations parallelized
        parfor k = 1:N_arcs
            sc = sim_config_arcs(k);
        
            % Propagate the orbit for a sigma value of zero
            ode_func = @(t, y) bcir4bp_stm(t, y, sc);
            
            sv_0 = [propagations(k).x0_k; reshape(eye(N), [N^2, 1])];
            ss_k = ode89(ode_func, [0, sim_config_arcs(k).prop_time], sv_0, opts);
        
            STM_aug_fk = reshape(ss_k.y(7:end, end), [N, N]);
        
            propagations(k).xf_k = ss_k.y(1:6, end);
            propagations(k).STM_fk = STM_aug_fk(1:6, 1:6);
        end

        % Loop over arcs to construct the constraint vector and DF matrix
        for k = 1:1:N_arcs
            index_range_a = 6*k-5:6*k;
            index_range_b = 6*k+1:6*k+6;
        
            if k < N_arcs
                % Enforce continuity from one node to the next
                F(index_range_a) = propagations(k).xf_k - X(index_range_b);
        
                % Constraint sensitivity with X
                DF(index_range_a, index_range_a) = propagations(k).STM_fk;
                DF(index_range_a, index_range_b) = -eye(6);
            else
                % Enforce continuity from the last node to the first
                F(index_range_a) = propagations(k).xf_k - X(1:6);
        
                DF(index_range_a, index_range_a) = propagations(k).STM_fk;
                DF(index_range_a, 1:6) = -eye(6);
            end
        end
        
        % Report progress on convergence to console
        fprintf("\t %d/%d, |F|:  %1.4e \n", q, q_max, norm(F))

        % Solve for change in X - important:  DF is treated as sparse.
        % Why?  I have no idea
        delta_x = sparse(DF) \ F;
    
        % Compute sigma stepsize adaptation factors
        if q == 1 % If this is the first rev through the N-R process
            first_SL = norm(delta_x); % How big was our first step?
        end
    
        if q == 2 % If this is the second rev through the N-R process
            contr_rate = norm(delta_x) / first_SL; % Stepsize changing fast?
        end    
    
        % If F is too big or small, break
        if norm(F) < convergence_tolerance
            break
        elseif norm(F) > 0.5
            q = q_max;  % Pretend that the loop has reached max iterations
            break
        end
        
        % If 'F' is neither converged nor diverged, apply the change in X 
        X = X - delta_x;
    end
    % ----- Newton-Raphson ----- %

    %% Plot intermediate orbits
    if and(plot_converged_orbits, mod(M, plot_converged_orbits) == 0)
        for k = 1:1:N_arcs
            sim_config_arcs(k).sigma = last_sigma;
            propagations(k).x0_k = last_X(6*k-5:6*k);
        end

        % Plot the original orbit for later comparison with PAC results
        figure(1); hold on; axis equal; grid on;
        
        this_orbit_color = rand([3, 1]);

        % Run integrations parallelized
        for k = 1:N_arcs
            sc = sim_config_arcs(k);
        
            % Propagate the orbit for a sigma value of zero
            ode_func = @(t, y) bcir4bp_stm(t, y, sc);
            
            sv_0 = [propagations(k).x0_k; reshape(eye(N), [N^2, 1])];
            ss_k = ode89(ode_func, [0, sim_config_arcs(k).prop_time], sv_0, opts);
        
            plot3(ss_k.y(1, :), ss_k.y(2, :), ss_k.y(3, :), 'Color', this_orbit_color)
        end

        script_pause = true;
    end

    %% Continue on with steplength adaptation and taking a step in sigma
    % Nominal values for steplength adaptation, compared against actuals
    nominal_contr_rate = 1.5; % Larger => bigger ds
    nominal_first_SL = 0.08;

    % Add some hard-coded bounds to cap/lower bound the adaptation factor
    maximum_adaptation = 2;
    minimum_adaptation = 0.5;

    % Add some hard-coded bounds to keep the adaptive scheme from moving fast/slow
    minimum_steplength = min_ds;
    maximum_steplength = 0.1;

    % If the chosen member was converged on iter #1, we need to set a
    % dummy contraction rate for that iteration.  Choose one that allows
    % for the other convergence metrics to be applied instead.
    if ~exist("contr_rate", "var")
        contr_rate = first_SL / nominal_first_SL * nominal_contr_rate;
    end

    % Two factors can be used to adapt our steplength
    adapt_metrics = [sqrt(contr_rate/nominal_contr_rate), ...
                     sqrt(first_SL/nominal_first_SL)];

    % Which of the two rooted ratios above is bigger?  Use it to limit ds
    largest_factor = max(adapt_metrics);

    % The more limiting of the two is chosen as our actual adaptation
    adapt_factor = max( min(largest_factor, maximum_adaptation), ...
                        minimum_adaptation);

    % A check to keep us from exceeding sigma = 1
    remaining_sigma = abs(1 - sigma);

    % Adapt ds using our control, but also upper/lower bound it
    ds = min(max(ds/adapt_factor, minimum_steplength), maximum_steplength);

    % Did the last convergence process fail?  If so, generate a new ds, retry
    if q >= q_max  
        fail_counter = fail_counter + 1;

        if fail_counter >= 100
            % Add a condition to die if working stepsize isn't found reasonably fast
            fprintf("Script failed! \n")
            break
        end

        % Generate a new decadally-uniform stepsize
        ds_exponent = rand * (ds_exponent_max - ds_exponent_min) + ds_exponent_min;

        % Choose the new ds, or most of the distance to '1' - whichever's smaller
        ds = min(10^(ds_exponent), 0.99 * remaining_sigma);
        
        % the sigma = 1 convergence is held to a higher tol - reset the
        % tol to the typical
        convergence_tolerance = 1e-6;

        fprintf("Last convergence failed, trying ds = %1.3e \n", ds)

        sigma = last_sigma + ds;
        X = last_X;

        continue
    end

    % If 'ds' is low for too long, we want to try larger stepsizes to
    % see if they still work
    if and(ds < 100 * min_ds, rand < (1/30))
        % If ds is within two orders of magnitude of min, shake up the
        % convergence process every thirty convergences on average

        % Generate a new decadally-uniform stepsize
        ds_exponent = rand * (ds_exponent_max - ds_exponent_min) + ds_exponent_min;

        % Choose the new ds, or the distance to '1' - whichever's smaller
        ds = min(10^(ds_exponent), remaining_sigma);

        fprintf("Shakeup occured, trying ds = %1.3e \n", ds)
    end

    % If we've reached this point, convergence has occured - reset fail ctr
    fail_counter = 0;

    % Report success to console
    fprintf("M: %d | Sigma:  %1.6f | ds:  %1.3e \n", M, sigma, ds)

    if sigma >= 1
        % Script has succeeded, break out of loop
        break
    elseif ds == remaining_sigma
        % Hold the last convergence @ sigma = 1 to a higher standard
        convergence_tolerance = 1e-10;
    end

    % Save last converged answer for posterity
    last_X = X;
    last_sigma = sigma;

    % Apply 'ds' to sigma in a way that doesn't exceed sigma = 1
    sigma = sigma + min(ds, remaining_sigma);
end