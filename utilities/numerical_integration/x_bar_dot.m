function xb_dot = x_bar_dot(th_s, xb, mu, m_s, a_s)
    % This function, using the characteristic mu of the three-body system,
    % the Cartesian state, and the solar angle returns the time
    % rate-of-change of the Cartesian state.
    
    x = xb(1); y = xb(2); z = xb(3);

    % The sun angle in the synodic frame is analytically calculable
    cth = cos(th_s);  sth = sin(th_s);

    d = sqrt((x+mu)^2 + y^2 + z^2);  % Distance to Earth, [nd]
    r = sqrt((x-1+mu)^2 + y^2 + z^2);  % Distance to Moon, [nd]
    r43 = sqrt((x-a_s*cth)^2 + (y-a_s*sth)^2 + z^2);  % Distance to Sun, [nd]
    
    xb_dot(1) = xb(4);
    xb_dot(2) = xb(5);
    xb_dot(3) = xb(6);
    
    % Similar to CR3BP plus one or two added terms
    xb_dot(4) = 2*xb(5)+x-(1-mu)*(x+mu)/d^3-mu*(x-1+mu)/r^3-m_s*(x-a_s*cth)/r43^3-m_s*cth/a_s^2;
    xb_dot(5) = -2*xb(4)+y-(1-mu)*y/d^3-mu*y/r^3-m_s*(y-a_s*sth)/r43^3-m_s*sth/a_s^2;
    xb_dot(6) = -(1-mu)*z/d^3-mu*z/r^3-m_s*z/r43^3;
end