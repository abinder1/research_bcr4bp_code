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
opts = odeset("RelTol",3e-14,"AbsTol",3e-14,"MaxStep",0.03);

peri_func = @(t,y) periapsis_pure(t, y, mu, 1, true);
opts_peri = odeset("RelTol",3e-14,"AbsTol",3e-14,"MaxStep",0.03,"Events",peri_func);

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

%% ODE45 function testing

dv = 37;  % Choose delta-v dataset we're interested in and load
input_data = compose('saved data/generated/periapsis_maps/ICs/pert_states_%d.mat', dv);
load(input_data{1});  % Loading in 'pert_states' structure

% Compose +dV and -dV into a single vector 
IC_set = [pert_states.state_neg pert_states.state_pos];

% In theory, run simulations at 100 distinct epochs denoted by starting Sun angle
N_sun_discr = 100;

% We're breaking up the N = 32000 fixed points per delta-v into chunks of 4000 sims (memory concerns)
sim_res(4000).q = [];
sim_res(4000).perix = [];
sim_res(4000).periy = [];

% Where do we want to store simulation results?
data_folder = compose('saved data/generated/periapsis_maps/result_chunks/%d/', dv);

k = 1;  % Index for starting Sun angle

for th_S0 = linspace(0, 2*pi, N_sun_discr)
    
    % Change the epoch by redefining the function handle
    ode_func = @(t,y)state_vec_derivs(t, y, mu_nd = mu, ...
                                            th_S0 = th_S0, ...
                                            m_s = m_s, ...
                                            a_s = a_s      );

    for chunk_num = 1:8  % = 32000 / 4000 = 8 chunks of sims
        chunk.idx_range = [4000*(chunk_num-1) + 1, 4000*chunk_num];
        chunk.th_S0 = th_S0;
        chunk.threshold = 0.15;
        
        % Get this chunk's subset of ICs
        chunk_ICs = IC_set(:, chunk.idx_range(1):chunk.idx_range(2));
        
        % We want to save periapses closer than 0.15 distance units from Earth
        % Prevents us from saving periapses too far from Earth.
        thr = chunk.threshold;

        parfor q = 1:10  % Parallelization for one chunk
            % -84 nondimensional time units is about 1 year into the past
            sol_struct = ode89(ode_func, [0, -84], chunk_ICs(:, q), opts_peri);
        
            % Keep finding periapses until one year has passed
            while sol_struct.x(end) > -84
                sol_struct = odextend(sol_struct, [], -84);
            end
    
            % After the above loop, all periapses are in sol_struct.ye
            % Which of this sim's periapses are below our chosen threshold?
            low_peris = vecnorm(sol_struct.ye(1:3, :) + [mu; 0; 0], 2) <= thr;
    
            % Save this sim's results to structs 'xe' and 'ye',
            % but only save results below our threshold
            sim_res(q).q = q + 4000*(chunk_num-1);
            sim_res(q).perix = sol_struct.xe(low_peris)
            sim_res(q).periy = sol_struct.ye(:, low_peris)
        end

        has_peris = false([1, 4000]);

        % ------------------------
        % This is a data filter
        for q = 1:1:length(has_peris)
            if ~isempty(sim_res(q).perix)
                has_peris(q) = true;
            end
        end

        % Filter out simulation entries that have no good periapses
        chunk.sim_res = sim_res(has_peris);
        % ------------------------

        % Save remaining chunk sims to file, indexed by Sun angle index and chunk #
        filename = compose('%d_%d_%d.mat', k, N_sun_discr, chunk_num);
        filename = fullfile(data_folder{1}, filename{1});

        save(filename, 'chunk')
    end

    % Increment Sun angle index by hand
    k = k + 1;
end