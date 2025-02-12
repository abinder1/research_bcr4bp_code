function sv_dot = bcir4bp_stm(delta_tau, sv, sim_config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the Bicircular Inclined Restricted Four-Body Problem (BCIR4BP),
% this function computes state derivatives and the state transition matrix
% in a format compatible with MATLAB's 'ode' suite.  Alongside the
% traditional STM, the partial derivatives of how the final state
% changes with the following four model parameters are also propagated:
%   1) The Sun's strength, \sigma
%   2) A change in the Earth/Moon - Sun semimajor axis, \tilde{a}_S
%
% Author:  Andrew Binder (2024)
%
% Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Handle both state-only and state and STM propagation
    if length(sv) == 6
        state = sv(1:6);
    else
        state = sv(1:6);
        STM = reshape(sv(7:end), [sim_config.N, sim_config.N]);
    end

    % Decompose quantities that come from the simulation config
    tau_e = sim_config.tau_e;
    mu = sim_config.mu;
    mutil_S = sim_config.mutil_S;
    atil_S = sim_config.atil_S;
    sigma = sim_config.sigma;
    B0 = sim_config.B0;
    TH0 = sim_config.TH0;
    INC = sim_config.INC;

    % Construct the angles at this particular \Delta \tau
    abs_time = delta_tau + tau_e;

    b_dot = sqrt((mutil_S + 1) / atil_S^3);

    theta = mod(TH0 + abs_time, 2*pi);
    B = mod(B0 - b_dot * abs_time, 2*pi);
    
    % Construct the C_{31B} direction cosine matrix
    cT = cos(theta);  sT = sin(theta);
    cI = cos(INC);  sI = sin(INC);
    cB = cos(B);  sB = sin(B);

    C3T = [cT, sT, 0; -sT, cT, 0; 0, 0, 1];
    C1I = [1, 0, 0; 0, cI, sI; 0, -sI, cI];
    C3B = [cB, sB, 0; -sB, cB, 0; 0, 0, 1];

    C_31B = C3T * C1I * C3B;

    % Earth -> satellite vector and unit vector
    rho_E = state(1:3) + [mu; 0; 0];
    rhohat_E = rho_E / norm(rho_E);

    % Moon -> satellite vector and unit vector
    rho_M = state(1:3) - [1 - mu; 0; 0];
    rhohat_M = rho_M / norm(rho_M);

    % Sun -> satellite vector and unit vector
    rho_S = state(1:3) + C_31B * [atil_S; 0; 0];  % Hardcoded vector is \tilde{d}
    rhohat_S = rho_S / norm(rho_S);

    % Acceleration of the Earth-Moon barycenter (EMBC)
    atil_EM = [mutil_S / atil_S^2; 0; 0];

    % Earth/Moon/Sun accelerations on the satellite wrt EMBC
    Atil_E = -(1 - mu) * rhohat_E / norm(rho_E)^2;
    Atil_M = -(mu) * rhohat_M / norm(rho_M)^2;
    Atil_S = -(mutil_S) * rhohat_S / norm(rho_S)^2;

    % Skew symmetric [\mathbbm{1}_3]_\times matrix from rotating frame
    skew_13 = [0, -1, 0; 1, 0, 0; 0, 0, 0];

    acceleration =  Atil_E ...                                          % Earth term
                  + Atil_M ...                                          % Moon term
                  + sigma * (Atil_S + C_31B * atil_EM) ...              % Sun-related terms
                  - 2 * skew_13 * state(4:6) - skew_13^2 * state(1:3);  % Terms from rotating frame

    sv_dot = [state(4:6); acceleration];

    if length(sv) > 6  % If we choose to integrate an STM
        A = zeros(7);

        A(1:3, 4:6) = eye(3);
        A(4:6, 4:6) = -2 * skew_13;

        % Jacobians of Earth/Moon/Sun acceleration terms
        dAtil_E_drho = (1 - mu) * (3 * (rhohat_E * rhohat_E') - eye(3)) / norm(rho_E)^3;
        dAtil_M_drho = mu * (3 * (rhohat_M * rhohat_M') - eye(3)) / norm(rho_M)^3;
        dAtil_S_drho = mutil_S * (3 * (rhohat_S * rhohat_S') - eye(3)) / norm(rho_S)^3;

        % Total partial w.r.t. spacecraft position
        A(4:6, 1:3) = dAtil_E_drho + dAtil_M_drho + sigma * dAtil_S_drho - skew_13^2;

        % Partial derivative with respect to sigma
        A(4:6, 7) = Atil_S + C_31B * atil_EM;

        STM_dot = A * STM;

        sv_dot = [sv_dot; reshape(STM_dot, [sim_config.N^2, 1])];
    end
end
