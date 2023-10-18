function [sv_dot] = state_vec_derivs(t, sv, nv_args)
    % This function serves as an unpackaging and packaging func.
    % for the tau-derivatives of both the Cartesian state in the
    % BCR4BP. Unpackage the row vector containing the state, find
    % each component's respective tau-derivative, and repackage
    % for ODE45.

    %  Notation:  Scheuerle's 2021 Thesis

    arguments % Allow for name-value argument definitions
        t double
        sv double
        nv_args.th_S0 double
        nv_args.mu_nd double
        nv_args.m_s double
        nv_args.a_s double
    end

    % Unpacking and naming for clarity
    x_bar = sv(1:6);

    % The sun angle in the synodic frame is analytically calculable
    th_s = nv_args.th_S0 + t * ( sqrt( (nv_args.m_s + 1) / nv_args.a_s^3 ) - 1 );
    
    xb_dot = x_bar_dot(th_s, x_bar, nv_args.mu_nd, ...
                                    nv_args.m_s, ...
                                    nv_args.a_s         );

    % Repackaged into a column vector
    sv_dot = xb_dot';
end