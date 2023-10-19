function Upsilon = BCR4BP_EM_PsuPot(xb, t, mu, m_s, a_s, th_S0)
    % This function computes the CR3BP pseudo-potential's value
    % at a point (x, y, z) ('a'-vector position).

    th_s = th_S0 + t * ( sqrt( (m_s + 1) / a_s^3 ) - 1 );

    % The sun angle in the synodic frame is analytically calculable
    cth = cos(th_s);  sth = sin(th_s);

    x = xb(1); y = xb(2); z = xb(3);

    d = sqrt((x+mu)^2 + y^2 + z^2);
    r = sqrt((x-1+mu)^2 + y^2 + z^2);
    r43 = sqrt((x-a_s*cth)^2 + (y-a_s*sth)^2 + z^2);  % Distance to Sun, [nd]

    Upsilon = (x^2 + y^2)/2 + (1-mu)/d + mu/r + m_s/r43 - m_s/a_s^2*(x*cth + y*sth);
end