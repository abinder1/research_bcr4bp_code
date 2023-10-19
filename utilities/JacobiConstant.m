function JC = JacobiConstant(xb, mu)
    % This function takes in an 'a'-vector, 'xb' and CR3BP nd \mu
    % value and returns the Jacobi constant of the state
    % represented by 'a'.
    
    JC = 2 * CR3BP_PsuPot(xb, mu) - dot(xb(4:6), xb(4:6));
end