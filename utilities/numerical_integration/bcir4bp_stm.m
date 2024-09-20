function sv_dot = bcir4bp_stm(tau, sv, sigma, M0, inc, RAAN, mu, ae, mu_S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the Bicircular Inclined Restricted Four-Body Problem (BCIR4BP),
% this function computes state derivatives in a format compatible with
% MATLAB's 'ode' suite.
%
% Author:  Andrew Binder (2024)
%
% Inputs:
%   tau = [T](1x1)<float> | The time elapsed since the simulation epoch
%   sv = [L, L/T](6x1)<float> | The nondimensionalized position and
%       velocity of a spacecraft flying within the model, expressed in
%       the CR3BP-typical rotating frame
%   sigma = [](1x1)<float> | A system configuration scalar that can
%       tune the effects caused by the Sun's gravity acting on the
%       model.  When sigma = 0, the model is identical to the CR3BP.
%       When sigma = 1, the model is the BCIR4BP with the Sun acting at
%       full-strength.  When sigma \in (0, 1), Sun effects are at
%       partial strength.
%   M0 = [rad](1x1)<float> | The argument of latitude of the Moon in its
%       circular orbit about the Earth, as measured at epoch and against
%       the Earth-Moon-Sun barycenter-centric inertial frame
%   inc = [rad](1x1)<float> | The constant inclination of the Moon's circular
%       orbit about the Earth, as measured against the Earth-Moon-Sun
%       barycenter-centric inertial frame
%   RAAN = [rad](1x1)<float> | The constant inclination of the Moon's circular
%       orbit about the Earth, as measured against the Earth-Moon-Sun
%       barycenter-centric inertial frame
%   mu = [](1x1)<float> | The dimensional standard gravitational
%       parameter (SGP) of the Moon, divided by the sum of the dimensional 
%       values of the Earth and the Moon's SGP.
%   ae = [L](1x1)<float> | The nondimensionalized semi-major axis of the
%       Earth's orbit about the Sun.
%   mu_S = [L^3 / T^2](1x1)<float> | The nondimensionalized SGP of the Sun.
%
% Outputs:
%   sv_dot = [L/T, L/T^2](6x1)<float> | The derivatives of each state
%       quantity with respect to tau.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-allocate for both speed, and to ensure that this vector has
    % the right shape
    sv_dot = zeros([6, 1]);

    % State vector unpacking
    satellite_position = sv(1:3);

    % This term equals 1 - \gamma, when nondimensionalized
    OmG = mu_S / (mu_S + 1);

    % Get the Moon's AoL and the RAAN - Earth MA diff. at the current tau
    M = M0 + tau;
    B = RAAN - sqrt(mu_S/ae^3) * tau;

    % Construct appropriate simple-rotation direction cosine matrices...
    C3M = [cos(M), sin(M), 0; -sin(M), cos(M), 0; 0, 0, 1];
    C1I = [1, 0, 0; 0, cos(inc), sin(inc); 0, -sin(inc), cos(inc)];
    C3B = [cos(B), sin(B), 0; -sin(B), cos(B), 0; 0, 0, 1];

    % ... and DCM partial derivatives
    K2B = [-cos(B), -sin(B), 0; sin(B), -cos(B), 0; 0, 0, 0];

    % Construct compound DCMs, DCM partial derivatives
    CB = C3M * C1I * C3B;
    Kalfa = C3M * C1I * K2B;

    % The Earth-Moon barycenter's distance to the Sun, in R_EM coords
    eps_vec = [-ae; 0; 0];

    % Get spacecraft positions relative to:
    %   1) The Sun, written in R_EM coordinates
    %   2) The Earth, written in M_EM coordinates
    %   3) The Moon, written in M_EM coordinates
    DeltaS = CB' * satellite_position - eps_vec;
    DeltaE = satellite_position - [-mu; 0; 0];
    DeltaM = satellite_position - [1 - mu; 0; 0];

    % Vector norms, squared
    DeltaS2 = dot(DeltaS, DeltaS);
    DeltaE2 = dot(DeltaE, DeltaE);
    DeltaM2 = dot(DeltaM, DeltaM);

    % Gravity forces, expressed in frames listed above
    A_S = -mu_S * DeltaS / DeltaS2^(3 / 2);
    A_E = -(1 - mu) * DeltaE / DeltaE2^(3 / 2);
    A_M = -mu * DeltaM / DeltaM2^(3 / 2);

    % Kinematical contribs. from the Moon's orbit, in M_EM coordinates
    moon_orbit_kinematic_contrib = [sv(1) + 2 * sv(5); sv(2) - 2 * sv(4); 0];

    % Kinematical contribs. from the Earth's orbit, in R_EM coordinates
    % NOTE: This term is the negative of the term from documentation,
    %       and also divided out a common term of 'ae'
    earth_orbit_kinematical_contrib = mu_S * (OmG / ae^2) * Kalfa * [1; 0; 0];

    % Synodic acceleration, expressed in M_EM coordinates
    xdd_MeM = sigma * (CB * A_S - earth_orbit_kinematical_contrib) ...
            + A_E + A_M + moon_orbit_kinematic_contrib;

    % Packaging results for MATLAB
    sv_dot(1:3) = sv(4:6);
    sv_dot(4:6) = xdd_MeM;
end
