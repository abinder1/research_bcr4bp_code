function [f, terminate, dir] = periapsis_pure(~, y, mu, prim_num, stop)
    % This more complicated event function triggers on one of
    % three conditions:
    %
    %       1) A periapsis is reached
    %       2) P1's x-value is reached
    %       3) P2's x-value is reached
    %
    % Only condition (1) can serve as stopping and only when
    % selected by the user.  dir = [-1 -1 -1] catches only
    % periapses, dir = [1 1 1] would catch apopases and [0 0 0]
    % would catch both.

    if prim_num == 1
        f = dot([y(1) + mu, y(2), y(3)], [y(4), y(5), y(6)]);
    elseif prim_num == 2
        f = dot([y(1) - 1 + mu, y(2), y(3)], [y(4), y(5), y(6)]);
    end

    if stop == true
        terminate = 1;
    else
        terminate = 0;
    end
    
    dir = -1;
end

