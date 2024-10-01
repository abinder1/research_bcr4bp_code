function sv_dot = bcir4bp_angles_stm(tau, sv, nv_args)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the Bicircular Inclined Restricted Four-Body Problem (BCIR4BP),
% this function computes state derivatives and the state transition matrix
% in a format compatible with MATLAB's 'ode' suite.
%
% Author:  Andrew Binder (2024)
%
% Inputs:
%   tau = [T](1x1)<float> | The time elapsed since the simulation epoch
%   sv = [L, L/T](6x1)<float> | (IFF stmenabled == false) The
%       nondimensionalized position and velocity of a spacecraft flying within
%       the model, expressed in the CR3BP-typical rotating frame.
%   sv = [L, L/T](42x1)<float> | (IFF stmenabled == true) The
%       nondimensionalized position and velocity of a spacecraft flying within
%       the model, expressed in the CR3BP-typical rotating frame, with a
%       linearly-indexed copy of the state transition matrix appended to the
%       end.
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
%   stmenabled = [](1x1)<boolean> | A flag that instructs the function
%       to also propagate the BCIR4BP's STM (and associated
%       configuration partials, used in homotopies)
%
% Outputs:
%   sv_dot = [L/T, L/T^2](6x1)<float> | The derivatives of each state
%       quantity with respect to tau.
%   sv_dot = [L/T, L/T^2](42x1)<float> | The derivatives of each state
%       quantity with respect to tau, with the tau-derivative of the STM
%       (a matrix) linearly-indexed and appended to the end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allow for name-value argument definitions for code clarity
    arguments
        tau double
        sv double
        nv_args.earth_moon_massparam double
        nv_args.earth_sma_nondim double
        nv_args.sun_sgp_nondim double
        nv_args.sun_effect_slider double = 1.0
        nv_args.moon_arglat_at_epoch double = 0.0
        nv_args.moon_inclination double = 0.898; % Moons incl. is 5.145 degrees
        nv_args.moon_right_ascension double = 0.0
        nv_args.stm_enabled logical = true
    end

    % Unpack name-value arguments into more usable variables
    mu = nv_args.earth_moon_massparam;
    ae = nv_args.earth_sma_nondim;
    mu_S = nv_args.sun_sgp_nondim;
    sigma = nv_args.sun_effect_slider;
    M0 = nv_args.moon_arglat_at_epoch;
    inc = nv_args.moon_inclination;
    RAAN = nv_args.moon_right_ascension;
    stmenabled = nv_args.stm_enabled;

    % Pre-allocate for both speed, and to ensure that this vector has
    % the right shape
    if stmenabled == false
        sv_dot = zeros([6, 1]);
    else
        sv_dot = zeros([106, 1]);

        STM = reshape(sv(7:106), [10, 10]);
    end

    % State vector unpacking
    satellite_position = sv(1:3);

    % This term equals 1 - \gamma, when nondimensionalized
    OmG = mu_S / (mu_S + 1);

    % Get the Moon's AoL and the RAAN - Earth MA diff. at the current tau
    M = M0 + tau;
    B = RAAN - sqrt((mu_S + 1)/ae^3) * tau;

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

    if stmenabled == true
        Amat = zeros(10);

        K3 = [-1, 0, 0; 0, -1, 0; 0, 0, 0];
        K4 = [0, 1, 0; -1, 0, 0; 0, 0, 0];

        sun_tensor = mu_S * (3 * (DeltaS * DeltaS') - DeltaS2 * eye(3));
        earth_tensor = (1 - mu) * (3 * (DeltaE * DeltaE') - DeltaE2 * eye(3));
        moon_tensor = mu * (3 * (DeltaM * DeltaM') - DeltaM2 * eye(3));

        % The Jacobians of each force term with respect to the s/c position
        dAS_dpos = sigma * CB * sun_tensor * CB' * DeltaS2^(-5 / 2);
        dAE_dpos = earth_tensor * DeltaE2^(-5 / 2);
        dAM_dpos = moon_tensor * DeltaM2^(-5 / 2);

        % Partial of dynamics with respect to parameter sigma
        dAS_dsigma = CB * A_S - earth_orbit_kinematical_contrib;

        % Similarity transform of dAS_dpos
        G = CB * dAS_dpos * CB';

        K1M = [-sin(M), cos(M), 0; -cos(M), -sin(M), 0; 0, 0, 0];
        K1B = [-sin(B), cos(B), 0; -cos(B), -sin(B), 0; 0, 0, 0];
        K3B = [sin(B), -cos(B), 0; cos(B), sin(B), 0; 0, 0, 0];
        dC1I_dI = [0, 0, 0; 0, -sin(inc), cos(inc); 0, -cos(inc), -sin(inc)];

        dCB_dM0 = K1M * C1I * C3B;
        dCB_dOm = C3M * C1I * K1B;
        dCB_dIn = C3M * dC1I_dI * C3B;

        dKalfa_dM0 = K1M * C1I * K2B;
        dKalfa_dOm = C3M * C1I * K3B;
        dKalfa_dIn = C3M * dC1I_dI * K2B;

        % K3 and K4 matrices come from the kinematical contr. to dynamics
        Amat(1:3, 4:6) = eye(3);

        Amat(4:6, 1:3) = dAS_dpos + dAE_dpos + dAM_dpos - K3;
        Amat(4:6, 4:6) = 2 * K4;
        Amat(4:6, 7) = dAS_dsigma;

        Amat(4:6, 8) = sigma * (dCB_dM0 * A_S ...
                            + G * dCB_dM0' * satellite_position ...
                            - mu_S * (OmG / ae^2) * dKalfa_dM0 * [1; 0; 0]);

        Amat(4:6, 9) = sigma * (dCB_dOm * A_S ...
                            + G * dCB_dOm' * satellite_position ...
                            - mu_S * (OmG / ae^2) * dKalfa_dOm * [1; 0; 0]);

        Amat(4:6, 10) = sigma * (dCB_dIn * A_S ...
                            + G * dCB_dIn' * satellite_position ...
                            - mu_S * (OmG / ae^2) * dKalfa_dIn * [1; 0; 0]);

        STM_dot = Amat * STM;

        sv_dot(7:106) = reshape(STM_dot, [100, 1]);
    end
end
