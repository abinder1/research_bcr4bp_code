function U_star = CR3BP_PsuPot(xb, mu)
    % This function computes the CR3BP pseudo-potential's value
    % at a point (x, y, z) ('a'-vector position).

    x = xb(1); y = xb(2); z = xb(3);
    d = sqrt((x+mu)^2 + y^2 + z^2); r = sqrt((x-1+mu)^2 + y^2 + z^2);

    U_star = (x^2 + y^2)/2 + (1-mu)/d + mu/r;
end