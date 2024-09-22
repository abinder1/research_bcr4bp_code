% Code by Andrew Binder

%% MATLAB Initialization and MATLAB Constants Definition
clear; clc; close('all')

addpath(genpath('utilities'));  % Add folder and subfolders
addpath(genpath('saved data'));  % Add folder and subfolders

% Isometric view command angles
iso = [-30, asind(1/sqrt(3))];

% Nice manifold colors
coolblue = [0.5294, 0.8078, 0.9804];
solidblue = [0.2594, 0.4078, 0.9804];
coolred = [240,128,128] / 255;
solidred = [240,36,36] / 255;

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
opts = odeset("RelTol", 3e-14, "AbsTol", 3e-14, "MaxStep", 0.03);

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

%% BLT Propogation

% Load in one of the nice ballistic lunar transfers
load("bcr4bp\saved data\generated\periapsis_maps\selected_BLT_states\R21S2_BLT_2.mat")

% Properly choose 'epoch' Sun angle (RPO insertion @ t = 0)
th_S0_TOI = blt_prop_time * (sqrt((m_s + 1)/a_s^3) - 1);

% Prepare the MATLAB ODE suite with this calculated epoch
ode_func = @(t,y)state_vec_derivs(t, y, mu_nd = mu, ...
                                        th_S0 = th_S0_TOI, ...
                                        m_s = m_s, ...
                                        a_s = a_s           );

% Propogate these IC's from \tau = epoch -> \tau = epoch + sim_dur
sol_struct = ode89(ode_func, [0, -blt_prop_time], blt_IC, opts);

JC_vec = zeros([1, length(sol_struct.x)]);
HEM_vec = zeros([1, length(sol_struct.x)]);

for k = 1:1:length(sol_struct.x)
    JC_vec(k) = JacobiConstant(sol_struct.y(:, k), mu);
    HEM_vec(k) = Hamiltonian_EM(sol_struct.y(:, k), sol_struct.x(k), mu, m_s, a_s, th_S0_TOI);
end

figure(2); hold on; grid on; plot(sol_struct.x, HEM_vec)
xlabel('Time [nd]');  ylabel('Hamiltonian, Earth-Moon [nd]')

figure(3); hold on; grid on; plot(sol_struct.x, JC_vec)
xlabel('Time [nd]');  ylabel('Jacobi Constant [nd]')

figure(orbitviews.shalos);

% Plot this full simulation to visualize 'settling' period
settle = plot3(sol_struct.y(1, :), sol_struct.y(2, :), sol_struct.y(3, :), 'b');

% T = 15.8 for BLT1, 8.5 for BLT2
% Also produce a truncated sim from the Earth until lunar flyby
sol_struct_flyby = ode89(ode_func, [0, 8.5], blt_IC, opts);

% Plot traj. until lunar flyby + lunar close approach pt.
blt = plot3(sol_struct_flyby.y(1, :), sol_struct_flyby.y(2, :), sol_struct_flyby.y(3, :), 'g');
flyby = scatter3(sol_struct_flyby.y(1, end), sol_struct_flyby.y(2, end), sol_struct_flyby.y(3, end), 'filled', 'g');

% Also mark the RPO injection point
inject = scatter3(blt_FC(1), blt_FC(2), blt_FC(3), 'filled', 'ks');

figure(2); scatter(sol_struct_flyby.x(end), Hamiltonian_EM(sol_struct_flyby.y(:, end), sol_struct_flyby.x(end), mu, m_s, a_s, th_S0_TOI))
figure(3); scatter(sol_struct_flyby.x(end), JacobiConstant(sol_struct_flyby.y(:, end), mu))

figure(orbitviews.shalos);

% Simulate the 'RPO' (unconverged in BCR4BP) for vis. purposes
blt_FC(4:6) = (1 - 0.1 / norm(blt_FC(4:6) / v_star)) * blt_FC(4:6);
                
% Change the epoch by redefining the function handle
ode_func = @(t,y)state_vec_derivs(t, y, mu_nd = mu, ...
                                        th_S0 = 0, ...
                                        m_s = m_s, ...
                                        a_s = a_s      );

% And integrate the 'RPO'
sol_struct_rpo = ode89(ode_func, [0, 3.5 * pi], blt_FC, opts);

% Plot the 'RPO'
orbit = plot3(sol_struct_rpo.y(1, :), sol_struct_rpo.y(2, :), sol_struct_rpo.y(3, :), 'r');

%% Initialize a trajectory figure to animate within
% Orbital figure initialize - lib. points and E/M
orbitviews.inertial = figure();
axis equal; hold on; % grid on; hold on;

xlabel("x-distance from Earth [km]")
ylabel("y-distance from Earth [km]")
zlabel("z-distance from Earth [km]")
% End orbital figure initialize

orbitviews.inertial.Position = [140, 120, 2^9, 2^9]

% Define custom figure view for optimal viewing
% view([-125, 20]);  % Good for BLT 1
view([-175, 25]);  % Good for BLT 2

orbitviews.inertial.Children.Position = [0.1, 0, 0.85, 1]
orbitviews.inertial.Children.Color = 'none'
orbitviews.inertial.Children.XColor = 'none'
orbitviews.inertial.Children.YColor = 'none'
orbitviews.inertial.Children.ZColor = 'none'

%% Animation parameters

time_frame = [blt_prop_time * t_star / 86400, 0]; % Dimensional msd [epoch end_day]
days_ps = 3;  % Days per second
mv_fps = 60;  % Frames per second

time_frame_sec = time_frame * 86400;  % Translate the time frame into seconds
dim_sspf = days_ps * 86400 / mv_fps;  % Figure out how many simulated seconds there are per frame
num_frames = ceil((time_frame_sec(2) - time_frame_sec(1)) / dim_sspf)  % Time vector at appropriate spacing

anim_tvec = linspace(time_frame_sec(1), time_frame_sec(2), num_frames);
anim_tvec_nd = anim_tvec / t_star;

%% Pull animation data from the solution structure
nd_pos_rot = deval(sol_struct, anim_tvec / t_star - blt_prop_time, 1:3);
nd_vel_rot = deval(sol_struct, anim_tvec / t_star - blt_prop_time, 4:6);

% Extract and dimensionalize rotating state, time vectors
dim_state_rot = [l_star * nd_pos_rot; v_star * nd_vel_rot];
dim_state_inert = zeros(size(dim_state_rot));

% Get the state of the Moon
moon_state_rot = repmat([l_star * (1-mu); 0; 0; 0; 0; 0], [1, length(anim_tvec)]);
moon_state_inert = zeros(size(moon_state_rot));

% Get the state of the Moon
earth_state_rot = repmat([-l_star * mu; 0; 0; 0; 0; 0], [1, length(anim_tvec)]);
earth_state_inert = zeros(size(earth_state_rot));

% Get sun angle time history
th_s_history = anim_tvec_nd * ( sqrt( (m_s + 1) / a_s^3 ) - 1 );

% Get the state of the Moon
sun_vec_rot = repmat([0; 0; 0; 0; 0; 0], [1, length(anim_tvec)]);
sun_vec_rot(1, :) = 1.25 * l_star * cos(th_s_history);
sun_vec_rot(2, :) = 1.25 * l_star * sin(th_s_history);
sun_vec_inert = zeros(size(sun_vec_rot));

% Transform rotating state to dimensional ECI state
for k = 1:1:length(anim_tvec)
    dim_state_inert(:, k) = Ri_to_Ei(dim_state_rot(:, k), mu, 1/t_star, anim_tvec(k));
    moon_state_inert(:, k) = Ri_to_Ei(moon_state_rot(:, k), mu, 1/t_star, anim_tvec(k));
    earth_state_inert(:, k) = Ri_to_Ei(earth_state_rot(:, k), mu, 1/t_star, anim_tvec(k));
    sun_vec_inert(:, k) = Ri_to_Ei(sun_vec_rot(:, k), mu, 1/t_star, anim_tvec(k));
end

%% Begin animation with answers from propogation
tail_length = 20;  % How many frames do you want the tail to persist
long_tail_length = num_frames;

F(num_frames) = struct('cdata',[],'colormap',[]);
v = VideoWriter('animations/BLT_2_anim_blacktail.mp4', 'MPEG-4'); v.FrameRate = mv_fps;

open(v); max_oop = max(dim_state_inert(3, :))

bound_UL = max(max(dim_state_inert(1:3, :), [], 2), [5; 5; 0]*10^5);
bound_LL = min(min(dim_state_inert(1:3, :), [], 2), [-5; -5; 0]*10^5);

xlim([bound_LL(1), bound_UL(1)]);
ylim([bound_LL(2), bound_UL(2)]);
zlim([bound_LL(3), bound_UL(3)]);

% Good views:  [-95, 17] for BLT, [-122, -8] for settling

sunvec_scale = -0.16;

for k = 1:1:length(anim_tvec)
    tail_start = max(k, 0); tail_end = max(k - tail_length - 1, 1);
    tail_pts = tail_end:1:tail_start;

    ltail_start = max(k, 0); ltail_end = max(k - long_tail_length - 1, 1);
    ltail_pts = ltail_end:1:ltail_start;

    moon = moon_state_inert(1:3, k);
    this_pt = dim_state_inert(1:3, k);
    earth = earth_state_inert(1:3, k);

    sunvec_basepoint = sun_vec_inert(1:3, k);
    sunvec_dispscale = sunvec_scale * sun_vec_inert(1:3, k);

    this_pt_tail = dim_state_inert(1:3, tail_pts);
    this_pt_ltail = dim_state_inert(1:3, ltail_pts);
    moon_tail = moon_state_inert(1:3, tail_pts);
    
    title_txt = "Days Until Epoch: " + compose("%5.2f", -anim_tvec(k) / 86400);

    if k == 1 % Movie initialize
        inertial_plots.earth = surf(X_sphere * 6378 + earth(1), Y_sphere * 6378 + earth(2), Z_sphere * 6378 + earth(3), ...
                    'EdgeColor', '#0047AB', 'FaceColor', '#6495ED');
        inertial_plots.moon = surf(X_sphere * 1738 + moon(1), Y_sphere * 1738 + moon(2), Z_sphere * 1738 + moon(3), ...
            'EdgeColor', '#696969', 'FaceColor', '#848884');

        mt_plot = plot3(moon_tail(1,:), moon_tail(2,:), moon_tail(3,:), '--');

        satellite = scatter3(this_pt(1), this_pt(2), this_pt(3), 'ks');
        slt_plot = plot3(this_pt_ltail(1,:), this_pt_ltail(2,:), this_pt_ltail(3,:), 'k');
        st_plot = plot3(this_pt_tail(1,:), this_pt_tail(2,:), this_pt_tail(3,:), 'r');

        sundir = quiver3(sunvec_basepoint(1), sunvec_basepoint(2), sunvec_basepoint(3), ...
                         sunvec_dispscale(1), sunvec_dispscale(2), sunvec_dispscale(3), 0);

        set(sundir, {"Color", "LineWidth", "MaxHeadSize"}, {[1.00,0.41,0.16], 4, 4})

        title(title_txt)

        txt1 = 'Days/Second = ' + compose("%4.1f", days_ps);
        txt2 = 'FPS = ' + compose("%2d", mv_fps);

        subtitle(txt1 + ', ' + txt2)
    else
        set(inertial_plots.moon, {'XData', 'YData', 'ZData'}, ...
            {X_sphere * 1738 + moon(1), Y_sphere * 1738 + moon(2), Z_sphere * 1738 + moon(3)});
        set(inertial_plots.earth, {'XData', 'YData', 'ZData'}, ...
            {X_sphere * 6378 + earth(1), Y_sphere * 6378 + earth(2), Z_sphere * 6378 + earth(3)});
        set(sundir, {'XData', 'YData', 'UData', 'VData'}, ...
            {sunvec_basepoint(1), sunvec_basepoint(2), sunvec_dispscale(1), sunvec_dispscale(2)});

        set(mt_plot, {'XData', 'YData', 'ZData'}, {moon_tail(1,:), moon_tail(2,:), moon_tail(3,:)});
        set(satellite, {'XData', 'YData', 'ZData'}, {this_pt(1), this_pt(2), this_pt(3)});
        set(slt_plot, {'XData', 'YData', 'ZData'}, {this_pt_ltail(1,:), this_pt_ltail(2,:), this_pt_ltail(3,:)});
        set(st_plot, {'XData', 'YData', 'ZData'}, {this_pt_tail(1,:), this_pt_tail(2,:), this_pt_tail(3,:)});

        orbitviews.inertial.Children.Title.String = title_txt;
    end

    F(k) = getframe(orbitviews.inertial);

    writeVideo(v, F(k));
end
close(v);

fig = figure;
movie(fig, F, 1, mv_fps)