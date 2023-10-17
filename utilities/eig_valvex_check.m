function largest_err = eig_valvex_check(monodromy, V, D)
    % Check to see that eigenvalues and eigenvectors are paired properly
    % (lol they are that's how I designed the code)
    eig_valvex_check = zeros(1, 6);
    
    for m = 1:1:6
        eig_valvex_check(m) = norm(mtimes(monodromy, V(:, m)) - D(m) * V(:, m));
    end
    
    largest_err = max(eig_valvex_check);
end

