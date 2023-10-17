function [sv_dot] = conic_propogator(~, sv, mu_dim)
    % This function takes some 'A'-state and STM (packaged as a
    % 42x1 vector 'sv') and returns the derivatives in a form
    % appropriate for MATLAB's 'ode' suite.  Function intended
    % for dimensional states and \mu.

    % Unpackaging of 'A'-state
    rv = sv(1:3); vv = sv(4:6); r = norm(rv);
    x = rv(1); y = rv(2); z = rv(3);
    
    sv_dot = zeros(42, 1);  % Preallocate space for sv's deriv.

    % Compute derivatives based on EOMs
    sv_dot(1:3) = vv;
    sv_dot(4:6) = -mu_dim * rv / norm(rv)^3;

    % Unpacking of state transition matrix
    phi = reshape(sv(7:end)', [6,6]);

    % Preallocate submatrix of block matrix "A"
    A21 = zeros(3,3);

    % Compute elements of submatrix
    A21(1) = 3*mu_dim*x^2 / r^5 - mu_dim / r^3;
    A21(2) = 3*mu_dim*x*y / r^5;
    A21(3) = 3*mu_dim*x*z / r^5;
    A21(4) = A21(2);
    A21(5) = 3*mu_dim*y^2 / r^5 - mu_dim / r^3;
    A21(6) = 3*mu_dim*y*z / r^5;
    A21(7) = A21(3);
    A21(8) = A21(6);
    A21(9) = 3*mu_dim*z^2 / r^5 - mu_dim / r^3;

    % Construct block matrix "A"
    A = [zeros(3), eye(3); A21, zeros(3)];

    % Compute time derivative of STM and add to end of 'sv_dot'
    phi_dot = mtimes(A, phi);
    sv_dot(7:42) = phi_dot;
end

