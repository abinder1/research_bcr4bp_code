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
sim_config.sigma = 1.0;
sim_config.B0 = 0.0;
sim_config.TH0 = 0.0;

% Propagate the orbit for a sigma value of zero
ode_func = @(t, y) bcir4bp_stm(t, y, sim_config);

dv_size = 25 / (1000 * v_star);  % Perturbation delta-vee size
prop_time = -40;

REMIN = 6378 / l_star;
RMMIN = 1738 / l_star;
RMAX = 5;

surf(X_sphere * REMIN - mu, Y_sphere * REMIN, Z_sphere * REMIN, ...
     'EdgeColor', '#0047AB', 'FaceColor', '#6495ED');

surf(X_sphere * RMMIN + 1 - mu, Y_sphere * RMMIN, Z_sphere * RMMIN, ...
     'EdgeColor', '#0047AB', 'FaceColor', '#C0C0C0');

crash_event = @(t, y) perigee_function(t, y, mu, 1.01 * REMIN, 1.01 * RMMIN, RMAX);
opts_crash = odeset("RelTol", 1e-9, "AbsTol", 1e-9, "Events", crash_event);

M_max = 300;

oo_time = 0;  % Choose an initial value for the on-orbit time

d_oo_min = 1e-10;
d_oo_max = 0.1;

d_oo_time = d_oo_min;

dxf_target_min = 0.05;
dxf_target_max = 0.2;

break_flag = false;
final_position_storage = zeros([3, M_max]);

earth_close_approach_struct(12 * M_max).pass_states = 0;
earth_close_approach_struct(12 * M_max).return_TOF = 0;
earth_close_approach_struct(12 * M_max).final_oo_state = 0;
earth_close_approach_struct(12 * M_max).final_oo_time = 0;
earth_close_approach_struct(12 * M_max).dv_size = dv_size;
earth_close_approach_struct(12 * M_max).pass_distance = 0;

moon_close_approach_struct = earth_close_approach_struct;

QE = 1;
QM = 1;

for M = 1:1:M_max
    min_moon_pass = 1;
    min_earth_pass = 1;

    fprintf("M = %d | On-Orbit Time = %.3f / %.3f | doo = %.2e\n", M, oo_time, orbit.TIP, d_oo_time)

    [IC_pure, dx0_dt0] = estimate_initial_state_change(oo_time, spo_ss, dv_size);
    
    % ----- Actually perturb the initial condition ----- %
    IC_burn = IC_pure;
    IC_burn(4:6) = IC_burn(4:6) * (1 + dv_size / norm(IC_pure(4:6)));  % Perturb the velocity
    
    perturbed_sol_struct = ode89(ode_func, [0, prop_time], [IC_burn; reshape(eye(6), [36, 1])], opts_crash);
    
    % Pull any events that happened
    earth_perigee_events = find(perturbed_sol_struct.ie == 4);

    for k = 1:1:length(earth_perigee_events)
        pass_distance = norm(perturbed_sol_struct.ye(1:3, earth_perigee_events(k)) + [mu; 0; 0]);

        min_earth_pass = min(min_earth_pass, pass_distance);

        if pass_distance < 0.1
            earth_close_approach_struct(QE).final_oo_state = IC_pure;
            earth_close_approach_struct(QE).final_oo_time = oo_time;
            earth_close_approach_struct(QE).dv_size = dv_size;
    
            earth_close_approach_struct(QE).pass_states = perturbed_sol_struct.ye(1:6, earth_perigee_events(k));
            earth_close_approach_struct(QE).return_TOF = -perturbed_sol_struct.xe(earth_perigee_events(k));
            earth_close_approach_struct(QE).pass_distance = pass_distance;
    
            QE = QE + 1;
        end
    end

    moon_perigee_events = find(perturbed_sol_struct.ie == 5);

    for k = 1:1:length(moon_perigee_events)
        pass_distance = norm(perturbed_sol_struct.ye(1:3, moon_perigee_events(k)) - [1 - mu; 0; 0]);

        min_moon_pass = min(min_moon_pass, pass_distance);

        if pass_distance < 0.1
            moon_close_approach_struct(QM).final_oo_state = IC_pure;
            moon_close_approach_struct(QM).final_oo_time = oo_time;
            moon_close_approach_struct(QM).dv_size = dv_size;
        
            moon_close_approach_struct(QM).pass_states = perturbed_sol_struct.ye(1:6, moon_perigee_events(k));
            moon_close_approach_struct(QM).return_TOF = -perturbed_sol_struct.xe(moon_perigee_events(k));
            moon_close_approach_struct(QM).pass_distance = pass_distance;
    
            QM = QM + 1;
        end
    end

    % Pull some quantities for the adaptive stepsize scheme
    final_state = perturbed_sol_struct.y(1:6, end);
    final_time = perturbed_sol_struct.x(end);

    earth_distance = min(norm(final_state(1:3) + [mu; 0; 0]), min_earth_pass);
    moon_distance = min(norm(final_state(1:3) - [1 - mu; 0; 0]), min_moon_pass);

    dxf_target = min([earth_distance; moon_distance]) / 10;
    dxf_target = max(min(dxf_target, dxf_target_max), dxf_target_min);

    final_position_storage(:, M) = final_state(1:3);

    STMF = reshape(perturbed_sol_struct.y(7:42, end), [6 6]);
    dxf_norm_dt0 = norm(STMF * dx0_dt0);
    
    if M > 1
        adapt_metric = (0.9 * norm(final_state - last_final_state) + 0.1 * dxf_norm_dt0 * d_oo_time) / dxf_target;

        adapt_factor = max(min(max(adapt_metric), 10), 0.8);
        
        d_oo_time = d_oo_time / adapt_factor;

        % This senses a large change in propagation time, i.e. a new crash
        % or escaping from crashes
        prop_time_condition = abs(final_time / last_final_time - 1);

        if prop_time_condition > 0.05 % We've just gotten past a crash event
            d_oo_time = d_oo_min;
        end
    end

    if break_flag == true
        break
    end
    
    if (oo_time + d_oo_time) > orbit.TIP
        d_oo_time = orbit.TIP - oo_time;

        break_flag = true;
    end

    d_oo_time = max(min(d_oo_time, d_oo_max), d_oo_min);

    % Set up for the next loop
    oo_time = oo_time + d_oo_time;

    last_final_state = final_state;
    last_final_time = final_time;
end

final_position_storage(:, M+1:end) = [];
earth_close_approach_struct(QE:end) = [];
moon_close_approach_struct(QM:end) = [];

scatter3(final_position_storage(1, :), final_position_storage(2, :), final_position_storage(3, :), 'k.')

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

    terminal = [1; 1; 1; 0; 0];

    direction = [-1; -1; 1; -1; -1];
end