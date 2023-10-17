function [V, D] = eigensort(S)
    % Inbuilt function
    [V, D] = eig(S, 'vector');
    
    D = sort(D, 'descend', 'ComparisonMethod', 'real');
    mdx_set = zeros(1,6);
    
    for m = 1:1:6
        compvec = mtimes(S, V(:, m)) ./ V(:, m);
        mindex = [0, 0];
    
        for k = 1:1:6
            if sum(round(compvec - D(k), 8) == 0) > mindex(2)
                mindex = [k sum(round(compvec - D(k), 8) == 0)];
            end
        end
    
        mdx_set(m) = mindex(1);
    end
    
    % The eigenvalues are sorted and the ith eigenvalue corresponds
    % to the i'th column of V
    [~, indx] = sort(mdx_set);
    V = V(:, indx);
end