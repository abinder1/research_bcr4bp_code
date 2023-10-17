function Ei = Ri_to_Ei(Ri, Pe, mm, TSE)
    % Ei - column 6-vector of dimensional state, written in ECI
    % Pe - unsigned distance from barycenter to Earth
    % mm - the mean motion of the two primaries = 1/t_star
    % TSE - the duration of time that has elapsed since epoch
    % Ri - column 6-vector of dimensional state, written in rot.

    Q = Q_matrix(mm, TSE);
    e = [-Pe; 0; 0; 0; 0; 0];

    Ei = Q * (Ri - e);
end

