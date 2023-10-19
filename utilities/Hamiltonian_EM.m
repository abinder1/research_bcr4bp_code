function HEM = Hamiltonian_EM(xb, t, mu, m_s, a_s, th_S0)
    % This function takes in an 'a'-vector, 'xb' and CR3BP nd \mu
    % value and returns the Jacobi constant of the state
    % represented by 'a'.
    
    HEM = 2 * BCR4BP_EM_PsuPot(xb, t, mu, m_s, a_s, th_S0) - ...
            dot(xb(4:6), xb(4:6));
end