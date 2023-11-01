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
        f1 = dot([y(1) + mu, y(2), y(3)], [y(4), y(5), y(6)]);
    elseif prim_num == 2
        f1 = dot([y(1) - 1 + mu, y(2), y(3)], [y(4), y(5), y(6)]);
    end

    f2 = norm([y(1) + mu, y(2), y(3)]) - 6400 / 384400;
    f3 = norm([y(1) - 1 + mu, y(2), y(3)]) - 1750 / 384400;

    f = [f1; f2; f3];

    if stop == true
        terminate = [1; 1; 1];
    else
        terminate = [0; 1; 1];
    end
    
    dir = [-1; 0; 0];
end

