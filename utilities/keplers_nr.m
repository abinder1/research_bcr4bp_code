function [E, k] = keplers_nr(e, M, N_iter, tol)
    % Getting this f(E) right to a tolerance of 2.6e-9 ensures
    % that the Moon's position is correct to within 1m accuracy.
    % Due to the implementation, oftentimes that position is
    % known to a degree better than this tolerance.
    %
    % The eccentricity of the Moon's orbit is 0.0549006.

    for k = 1:1:N_iter
        if k == 1
            E = M;
        else
            E = E + f / df;
        end

        f = M - E + e*sin(E);
        df = 1 - e*cos(E);

        if abs(f) <= tol
            break
        end
    end
end