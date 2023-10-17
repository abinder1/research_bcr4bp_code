function w = fornberg_weights(m, z, s)
    %% Description of function
    % Full disclosure, I don't understand this implementation.
    % The algorithm and implementation comes directly from:
    %
    %   Fornberg, B. "Calculation of Weights in Finite Difference
    %       Formulas." SIAM Review, Vol. 40, No. 3, 1998,
    %       pp. 685–691.
    %
    %% Input(s) to the function
    % This function takes three inputs:
    %   'm' [int >= 0] - the deg. of deriv. you wish to take
    %       - e.g. if m = 2, you will get the acceleration
    %   'z' [scalar] - the loc where derivatives must be accurate
    %       - I'm pretty sure z must be in the range spanned by s
    %   's' [scalar list, length n] - grid pts of choice
    %       - Grid points can be arbitrary, but must ascend
    %
    %% Output(s) from the function
    % This function returns one output:
    %   'w' [scalar matrix, 'm+1' rows x 'n' columns] - The set
    %       of weights returned by the algorithm
    %       - The j'th row of 'w' are the weights for evaluating 
    %         the (j-1)'th derivative.  This can be calculated
    %         cleanly as follows:
    %           
    %           A = (j-1)'th row of w
    %           B = f(s) for all values in s
    %
    %           d^kf/dx^k = dot(A, B)
    %
    %% Input validation

    if ~issorted(s) % If 's' is not in ascending order
        error('Error calculating Fornberg weights - s must be in ascending order')
    end

%     if or(z > max(s), z < min(s)) % If 'z' is out of bounds
%         error('Error calculating Fornberg weights - z must lie in the span of s')
%     end

    if m >= length(s)  % Too high of a derivative requested
        error('Error calculating Fornberg weights - you need m + 1 points in s to calcuate the m''th derivative')
    end

    %% Recursive calculation of weights

    w = zeros(m+1, length(s));
    c = [1, 0, 0, s(1) - z, 0];
    
    w(1, 1) = 1;
    
    for i = 1:(length(s) - 1)
        mn = min(i, m);
        c(2) = 1;
        c(5) = c(4);
        c(4) = s(i+1) - z;
    
        for j = 0:(i-1)
            c(3) = s(i+1) - s(j+1);
            c(2) = c(2) * c(3);
    
            if j == (i-1)
                for k = mn:-1:1
                    w(k+1, i+1) = c(1) * (k * w(k, i) - c(5) * w(k+1, i)) / c(2);
                end
    
                w(1, i+1) = -c(1) * c(5) * w(1, i) / c(2);
            end
    
            for k = mn:-1:1
                w(k+1, j+1) = (c(4) * w(k+1, j+1) - k * w(k, j+1)) / c(3);
            end
    
            w(1, j+1) = c(4) * w(1, j+1) / c(3);
        end
        c(1) = c(2);
    end
end