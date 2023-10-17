function Q = Q_matrix(mm, TSE)
    % TSE - the duration of time that has elapsed since epoch
    % mm - the mean motion of the two primaries = 1/t_star

    tau = mm * TSE;

    I_C_R = [ cos(tau), sin(tau), 0; ...
             -sin(tau), cos(tau), 0; ...
              0,        0,        1]     ;
    
    I_Cdot_R = [sin(tau), -cos(tau), 0; ...
                cos(tau),  sin(tau), 0; ...
                0,         0,        0] * -mm;

    Q = [I_C_R' zeros(3); I_Cdot_R' I_C_R'];
end

