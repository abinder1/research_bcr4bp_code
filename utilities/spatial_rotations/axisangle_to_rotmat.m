function R = axisangle_to_rotmat(theta, u_hat)
    % This function creates a rotation matrix R.  This rotation
    % matrix is equivalent to a proper rotation R by an angle
    % theta about a axis specified by unit direction u_hat.  This
    % matrix can be used as Rq = q' where q is a column -
    % performing this operation will rotate q by theta radians
    % right-hand-rule about u_hat.

    ct = cos(theta); st = sin(theta);
    ux = u_hat(1); uy = u_hat(2); uz = u_hat(3); 
    
    R = zeros(3);
    R(1) = ct + ux^2*(1-ct);
    R(2) = ux*uy*(1-ct) + uz*st;
    R(3) = ux*uz*(1-ct) - uy*st;
    R(4) = ux*uy*(1-ct) - uz*st;
    R(5) = ct + uy^2*(1-ct);
    R(6) = uy*uz*(1-ct) + ux*st;
    R(7) = ux*uz*(1-ct) + uy*st;
    R(8) = uy*uz*(1-ct) - ux*st;
    R(9) = ct + uz^2*(1-ct);
end

