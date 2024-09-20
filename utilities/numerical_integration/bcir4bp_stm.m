function sv_dot = bcir4bp_stm(tau, sv, sigma, M0, inc, Omega, mu, aE, Bdot)
    sv_dot = zeros([6, 1]);

    % Get the angles at the current time
    M = M0 + tau;
    B = Omega - Bdot * tau;

    C3M = [cos(M), sin(M), 0; -sin(M), cos(M), 0; 0, 0, 1];
    C1I = [1, 0, 0; 0, cos(inc), sin(inc); 0, -sin(inc), cos(inc)];
    C3B = [cos(B), sin(B), 0; -sin(B), cos(B), 0; 0, 0, 1];

    K2B = [-cos(B), -sin(B), 0; sin(B), -cos(B), 0; 0, 0, 0];

    CB = C3M * C1I * C3B;
    Kalfa = C3M * C1I * K2B;

    eps_vec = [-aE; 0; 0];
    pos = sv(1:3);

    DeltaS = CB' * pos - eps_vec;
    DeltaE = pos - [-mu; 0; 0];
    DeltaM = pos - [1-mu; 0; 0];

    DeltaS2 = dot(DeltaS, DeltaS);
    DeltaE2 = dot(DeltaE, DeltaE);
    DeltaM2 = dot(DeltaM, DeltaM);

    A_Sgam = -(1 - sigma)^(-1) * deltaS / DeltaS2^(3/2);
    A_E = -(1 - mu) * deltaE / DeltaE2^(3/2);
    A_M = -mu * deltaM / DeltaM2^(3/2);

    k_vec = [sv(1) + 2 * sv(5); sv(2) - 2 * sv(4); 0];

    xdd_MeM = sigma * (CB * A_Sgam + Bdot^2 * Kalfa * eps_vec) + A_E + A_M + k_vec;

    sv_dot(1:3) = sv(4:6);
    sv_dot(4:6) = xdd_MeM;
end

