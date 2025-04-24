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

% Sphere mesh points for plotting of the primaries
[X_sphere, Y_sphere, Z_sphere] = sphere(10);

aE_dim = 149598023; % Earth's SMA about the Sun

mu = mu_moon / (mu_earth + mu_moon); % mu value for the Earth-Moon system
mu_oom = round(log10(mu));
mu_SF_known = 8; % We know 8 Earth sigfigs, 9 Moon sigfigs

% mu rounded to appropriate significance
mu = round(mu*10^(mu_SF_known - mu_oom - 1)) / 10^(mu_SF_known - mu_oom - 1);

muS = mu_sun / (mu_earth + mu_moon); % Sun's n.d. SGP [L^3 / T^2]

l_star = 384400;  % Equal to the average SMA of the Moon about EM Bary
t_star = sqrt(l_star^3/(mu_moon + mu_earth)); % Div. by 86400 > time in days
v_star = l_star / t_star; % Characteristic speed, in km/s

aE = aE_dim / l_star; % Nondimensionalize the Earth-Sun distance

%% Defintion of ODE options structures
opts = odeset("RelTol", 1e-12, "AbsTol", 1e-12);

%% Build out the right initial X vector via rev stacking
load('..\saved data\generated\l4_short_period.mat')

% #84 is the closest to Brian's initial condition
orbit = l4_short_period(84);

sim_config.tau_e = 0;
sim_config.mu = mu;
sim_config.mutil_S = muS;
sim_config.atil_S = aE;
sim_config.sigma = 0.0;
sim_config.B0 = 0.0;
sim_config.TH0 = 0.0;
sim_config.INC = deg2rad(5.145);

sim_config.simulate_STM = true;

% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_stm(t, y, sim_config);

sv0 = [orbit.ic; reshape(eye(6), [36, 1])];

spo_ss = ode89(ode_func, [0, orbit.TIP], sv0, opts);

% Plot the original orbit for later comparison with PAC results
figure(1); hold on; axis equal; grid on;

plot3(spo_ss.y(1, :), spo_ss.y(2, :), spo_ss.y(3, :), 'k')

%% Try to adjust d_oo_time to retain sensitivity
prop_time = -40;
sim_config.sigma = 1.0;

b_dot = sqrt((muS + 1) / aE^3);

REMIN = 6378 / l_star;
RMMIN = 1738 / l_star;
RMAX = 5;

surf(X_sphere * REMIN - mu, Y_sphere * REMIN, Z_sphere * REMIN, ...
     'EdgeColor', '#0047AB', 'FaceColor', '#6495ED');

surf(X_sphere * RMMIN + 1 - mu, Y_sphere * RMMIN, Z_sphere * RMMIN, ...
     'EdgeColor', '#0047AB', 'FaceColor', '#C0C0C0');

crash_event = @(t, y) perigee_function(t, y, mu, 1.01 * REMIN, 1.01 * RMMIN, RMAX);
opts_crash = odeset("RelTol", 1e-9, "AbsTol", 1e-9, "Events", crash_event);

d_oo_min = 1e-5;
d_oo_max = 0.01;

M_max = ceil(orbit.TIP / d_oo_min);

dxf_target = 0.005;

% final_position_storage = zeros([3, M_max]);

Qmax = 10000;

B0_space = [0, 1, 2, 3] * pi/2;
TH0_space = B0_space;
dv_space = [-50:25:50] / (1000 * v_star);

configuration_space = zeros([length(B0_space) * length(TH0_space) * length(dv_space), 3]);

k4 = 1;

for k1 = 1:1:length(B0_space)
    for k2 = 1:1:length(TH0_space)
        for k3 = 1:1:length(dv_space)
            configuration_space(k4, :) = [B0_space(k1), TH0_space(k2), dv_space(k3)];

            k4 = k4 + 1;
        end
    end
end

for k4 = 10:1:length(configuration_space)
    earth_close_approach_struct(Qmax).pass_states = 0;
    earth_close_approach_struct(Qmax).return_TOF = 0;
    earth_close_approach_struct(Qmax).pass_distance = 0;
    earth_close_approach_struct(Qmax).CA_B = 0;
    earth_close_approach_struct(Qmax).CA_TH = 0;

    earth_close_approach_struct(Qmax).M = 0;

    earth_close_approach_struct(Qmax).final_oo_state = 0;
    earth_close_approach_struct(Qmax).final_oo_time = 0;

    earth_close_approach_struct(Qmax).dv_size = 0;
    earth_close_approach_struct(Qmax).final_B = 0;
    earth_close_approach_struct(Qmax).final_TH = 0;

    earth_close_approach_struct = orderfields(earth_close_approach_struct);

    % theta = mod(TH0 + abs_time, 2*pi);
    % B = mod(B0 - b_dot * abs_time, 2*pi);
    oo_time = 0;  % Choose an initial value for the on-orbit time
    break_flag = false;
    d_oo_time = d_oo_min;

    sim_config.B0 = configuration_space(k4, 1);
    sim_config.TH0 = configuration_space(k4, 2);
    
    % Propagate the orbit for a sigma value of zero
    ode_func = @(t, y) bcir4bp_stm(t, y, sim_config);
    
    dv_size = configuration_space(k4, 3);  % Perturbation delta-vee size

    QE = 1;
    QE_log = 1;
    
    tic;
    
    for M = 1:1:M_max
        [IC_pure, dx0_dt0] = estimate_initial_state_change(oo_time, spo_ss, dv_size);
        
        % ----- Actually perturb the initial condition ----- %
        IC_burn = IC_pure;
        IC_burn(4:6) = IC_burn(4:6) * (1 + dv_size / norm(IC_pure(4:6)));  % Perturb the velocity
        
        perturbed_sol_struct = ode89(ode_func, [0, prop_time], [IC_burn; reshape(eye(6), [36, 1])], opts_crash);
        
        sensitive_states = zeros([6, 3]);
        max_dxf_norm_dt0 = zeros([1, 3]);
    
        while and(perturbed_sol_struct.ie(end) > 3, abs(perturbed_sol_struct.x(end)) < abs(prop_time))
            if perturbed_sol_struct.ie(end) == 4
                pass_distance = norm(perturbed_sol_struct.ye(1:3, end) + [mu; 0; 0]);
    
                if pass_distance < 0.15
                    earth_close_approach_struct(QE).final_oo_state = IC_pure;
                    earth_close_approach_struct(QE).final_oo_time = oo_time;
                    earth_close_approach_struct(QE).dv_size = dv_size;
            
                    earth_close_approach_struct(QE).pass_states = perturbed_sol_struct.ye(1:6, end);
                    earth_close_approach_struct(QE).return_TOF = -perturbed_sol_struct.xe(end);
                    earth_close_approach_struct(QE).pass_distance = pass_distance;
    
                    earth_close_approach_struct(QE).M = M;

                    earth_close_approach_struct(QE).final_B = sim_config.B0;
                    earth_close_approach_struct(QE).final_TH = sim_config.TH0;
                    earth_close_approach_struct(QE).CA_B = mod(sim_config.B0 - b_dot * perturbed_sol_struct.xe(end), 2*pi);
                    earth_close_approach_struct(QE).CA_TH = mod(sim_config.TH0 + perturbed_sol_struct.xe(end), 2*pi);
                         
                    QE = QE + 1;
    
                    if QE > Qmax
                        QE = 1;
    
                        fprintf("\t Saving Earth pass data (QE_log = %d)... \n\n", QE_log)
    
                        filename = sprintf("./perigee_data/B_%02d_T_%02d_DV_%02d_earth_%d.mat", round(sim_config.B0 * 10), round(sim_config.TH0 * 10), round(dv_size * (1000 * v_star)), QE_log);
                        save(filename, 'earth_close_approach_struct')
    
                        QE_log = QE_log + 1;
                    end
                end
            end
    
            dxf_norm_dt0 = norm(reshape(perturbed_sol_struct.y(7:42, end), [6 6]) * dx0_dt0);
            
            if any(max_dxf_norm_dt0 < dxf_norm_dt0)
                max_dxf_norm_dt0_aug = [max_dxf_norm_dt0, dxf_norm_dt0];
                sensitive_states_aug = [sensitive_states, perturbed_sol_struct.y(1:6, end)];
    
                [max_dxf_norm_dt0, max_dxf_indices] = sort(max_dxf_norm_dt0_aug, 'descend');
    
                max_dxf_norm_dt0 = max_dxf_norm_dt0(1:3);
                sensitive_states = sensitive_states_aug(:, max_dxf_indices(1:3));
            end
    
            perturbed_sol_struct = odextend(perturbed_sol_struct, [], prop_time);
        end
        
        dxf_norm_dt0 = norm(reshape(perturbed_sol_struct.y(7:42, end), [6 6]) * dx0_dt0);
    
        final_earth_distance = norm(perturbed_sol_struct.ye(1:3, end) + [mu; 0; 0]);
        final_moon_distance = norm(perturbed_sol_struct.ye(1:3, end) - [1 - mu; 0; 0]);
    
        closest_primary = min([final_earth_distance, final_moon_distance]);
        
        if and(any(max_dxf_norm_dt0 < dxf_norm_dt0), closest_primary < 0.2)
            max_dxf_norm_dt0_aug = [max_dxf_norm_dt0, dxf_norm_dt0];
            sensitive_states_aug = [sensitive_states, perturbed_sol_struct.y(1:6, end)];
    
            [max_dxf_norm_dt0, max_dxf_indices] = sort(max_dxf_norm_dt0_aug, 'descend');
    
            max_dxf_norm_dt0 = max_dxf_norm_dt0(1:3);
            sensitive_states = sensitive_states_aug(:, max_dxf_indices(1:3));
        end
    
        % final_position_storage(:, M) = perturbed_sol_struct.y(1:3, end);
        
        if M > 1
            sensitive_state_norms = vecnorm(sensitive_states - last_sensitive_states, 2);
    
            adapt_metric = (1.00 * max(sensitive_state_norms) + 0.00 * max(max_dxf_norm_dt0) * d_oo_time) / dxf_target;
    
            adapt_factor = min(max(adapt_metric, 1/2), 10);
            
            d_oo_time = d_oo_time / adapt_factor;
        end
    
        d_oo_time = max(min(d_oo_time, d_oo_max), d_oo_min);
    
        if break_flag == true
            break
        end
        
        if (oo_time + d_oo_time) > orbit.TIP
            d_oo_time = orbit.TIP - oo_time;
    
            break_flag = true;
        end
    
        elapsed_time = toc;
    
        M_remaining_expectation = round((orbit.TIP - oo_time) / d_oo_time);
        M_per_second = M / elapsed_time;
    
        time_remaining_expectation = M_remaining_expectation / M_per_second;
    
        el_hours = floor(elapsed_time / 3600);
        el_minutes = floor(elapsed_time / 60 - 60 * el_hours);
        el_seconds = floor(elapsed_time - 3600 * el_hours - 60 * el_minutes);
        duration_string = sprintf("%02d:%02d:%02d", el_hours, el_minutes, el_seconds);
    
        re_hours = floor(time_remaining_expectation / 3600);
        re_minutes = floor(time_remaining_expectation / 60 - 60 * re_hours);
        re_seconds = floor(time_remaining_expectation - 3600 * re_hours - 60 * re_minutes);
        remaining_string = sprintf("%02d:%02d:%02d", re_hours, re_minutes, re_seconds);
    
        fprintf("k4 = %d / %d | M = %d | On-Orbit Time = %.3f / %.3f | doo = %.2e | Elapsed: %s | Remaining: %s\n", k4, length(configuration_space), M, oo_time, orbit.TIP, d_oo_time, duration_string, remaining_string)
    
        % Set up for the next loop
        oo_time = oo_time + d_oo_time;
    
        last_sensitive_states = sensitive_states;
    end
    
    % final_position_storage(:, M+1:end) = [];
    % scatter3(final_position_storage(1, 1:M), final_position_storage(2, 1:M), final_position_storage(3, 1:M), 'k.')
    
    earth_close_approach_struct(QE:end) = [];
    
    fprintf("\t Saving Earth pass data \n\n")
    
    filename = sprintf("./perigee_data/B_%02d_T_%02d_DV_%02d_earth_%d.mat", round(sim_config.B0 * 10), round(sim_config.TH0 * 10), round(dv_size * (1000 * v_star)), QE_log);
    save(filename, 'earth_close_approach_struct')
end


%% Local function definitions
function [IC_pure, dx0_dt0] = estimate_initial_state_change(oo_time, sol_struct, dv_size)
    ode_func = sol_struct.extdata.odefun;

    IC_pure = deval(sol_struct, oo_time, 1:6);  % Pick an IC on the orbit
    
    sv_dot_0 = ode_func(0, [IC_pure; reshape(eye(6), [36, 1])]);
    
    dp0_dt0 = sv_dot_0(1:3);
    
    dvelfac_dt0 = -dv_size * dot(IC_pure(4:6), sv_dot_0(4:6)) / norm(IC_pure(4:6))^3;

    dv0_dt0 = dvelfac_dt0 * IC_pure(4:6) + (1 + dv_size / norm(IC_pure)) * sv_dot_0(4:6);

    dx0_dt0 = [dp0_dt0; dv0_dt0];
end

function [f, terminal, direction] = perigee_function(~, y, mu, REMIN, RMMIN, RMAX)
    pos = y(1:3);
    vel = y(4:6);

    earth_crash_function = norm(pos + [mu; 0; 0]) - REMIN;
    moon_crash_function = norm(pos - [1 - mu; 0; 0]) - RMMIN;
    escape_function = norm(pos) - RMAX;
    earth_pass = dot(pos + [mu; 0; 0], vel);
    moon_pass = dot(pos - [1 - mu; 0; 0], vel);

    f = [earth_crash_function; ...
         moon_crash_function; ...
         escape_function; ...
         earth_pass; ...
         moon_pass];

    terminal = [1; 1; 1; 1; 1];

    direction = [];
end