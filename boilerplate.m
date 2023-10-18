% Code by Andrew Binder

%% MATLAB Initialization and MATLAB Constants Definition
clear; clc; close('all')

addpath(genpath('utilities'));  % Add folder and subfolders
addpath(genpath('saved data'));  % Add folder and subfolders

% Red-blue colormap for Moon distance coloring
cmap_lunar = [linspace(1, 0, 1000); zeros(1, 1000); linspace(0, 1, 1000)]';

% Isometric view command angles
iso = [-30, asind(1/sqrt(3))];

% Nice manifold colors
coolblue = [0.5294, 0.8078, 0.9804];
solidblue = [0.2594, 0.4078, 0.9804];
coolred = [240,128,128] / 255;
solidred = [240,36,36] / 255;

% Sphere mesh points for plotting of the primaries
[X_sphere, Y_sphere, Z_sphere] = sphere(10);

%% Permit usage of NASA NAIF's MICE Toolkit
addpath(fullfile(pwd, './mice/src/mice/'))
addpath(fullfile(pwd, './mice/lib/'))

%% Constants of the Problem
mu_moon = 4902.8005821478;  % Both in km^3 / s^2
mu_earth = 398600.4415;
mu_sun = 132712440018;

l_star = 384400;  % Also equal to the SMA of the Moon
t_star = sqrt(l_star^3 / (mu_moon+mu_earth)); % Divide by 86400 to get time in days
v_star = l_star / t_star;

a_s = 149.598 * 10^6 / l_star;  % E-M Barycenter's SMA about the Sun, [nd]

mu = mu_moon / (mu_earth + mu_moon);  % mu value for the Earth-Moon system
m_s = mu_sun / (mu_earth + mu_moon);  % Nondimensionalized solar mass

mu_oom = round(log10(mu)); mu_SF_known = 8;  % We know 8 Earth sigfigs, 9 Moon sigfigs

% mu rounded to appropriate significance
mu = round(mu*10^(mu_SF_known-mu_oom-1)) / 10^(mu_SF_known-mu_oom-1);

%% Definition of option structures
opts = odeset("RelTol",1e-10,"AbsTol",1e-12);

%% Ephemeris file loading and access

% Get good absolute paths to the datafiles
planet_datafile = fullfile(pwd, './saved data/NASA NAIF/DE405AllPlanets.bsp');
leapseconds_datafile = fullfile(pwd, './saved data/NASA NAIF/naif0011.tls');

% Load in the appropriate datafiles to MICE
cspice_furnsh( planet_datafile )  % Gives Earth, Moon positions/velocities
cspice_furnsh( leapseconds_datafile )  % Gives millisecond timing data

%% Load data from a custom compressed solution structure file
load("saved data\CR3BP_EMROs\p2q1_3D_linking_family.mat")

%% Organize solution structure data using custom structure
% Solution structure to organize MS data
int_results = struct('x_0k', {}, 'epoch', {}, 'int_time', {}, 'x_fk', {}, 'stm_fk', {});

%% Propogate an initial condition
% For this initial condition
sv_k = [0.9043, 0, 0.0411, 0, -0.7721, 0]';

% Change the epoch by redefining the function handle
ode_func = @(t,y)state_vec_derivs(t, y, mu_nd = mu, ...
                                        th_S0 = 0, ...
                                        m_s = m_s, ...
                                        a_s = a_s      );

% And integrate the 'RPO'
sol_struct = ode89(ode_func, [0, 3.5 * pi], sv_k, opts);

%% Code to normalize figures
f = figure(1);

view([0, 90]) % Define the perspective, [az., el.]

% Modify the figure's size explicitly
f.Units = 'inches';
f.Position(3:4) = [6, 4]; % 6 inches wide by 4 inches tall

% Modify the figure's axes position for centering
f.Children.Position(2) = 0.14; % Center the axes vertically

% Font changing done at 'Axes' object, use TNR for journals
f.Children.FontName = 'Times New Roman';
f.Children.FontSize = 12;

% Specify axes limits explicitly
f.Children.XLim = [-400000, 400000] / l_star;
f.Children.YLim = [-310000, 310000] / l_star;
f.Children.ZLim = [-310000, 310000] / l_star;

% Set tick marks explicitly
f.Children.XTick = (-300000:150000:300000) / l_star;
f.Children.YTick = (-300000:150000:300000) / l_star;
f.Children.ZTick = (-300000:150000:300000) / l_star;

% Neatly written tick labels
f.Children.XTickLabel = compose('%3.1f', -3:1.5:3);
f.Children.YTickLabel = compose('%3.1f', -3:1.5:3);
f.Children.ZTickLabel = compose('%3.1f', -3:1.5:3);

% Axis labels
xlabel('x-distance from barycenter [km]');
ylabel('y-distance from barycenter [km]');
zlabel('z-distance from barycenter [km]');