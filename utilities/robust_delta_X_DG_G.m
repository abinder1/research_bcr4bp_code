function [delta_x, step_residual] = robust_delta_X_DG_G(DG, G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This utility robustly solves a linear system G = DG * delta_x using a
% variety of methods:
%   1) GMRES, with conditioning applied to DG and G
%   2) GMRES, unconditioned
%   3) Solving via row reduction
%   4) Solving using MATLAB inbuilt 'DG \ G'
%
% A residual for each answer is computed by taking it's solution dx_i,
% generating it's corresponding G vector 'g_i', and finding the solution
% with the smallest residual = norm(g_i - G)
%
% Author:  Andrew Binder (2024)
%
% Inputs:
%   G = [](nx1)<float> | The left-hand side of the above system
%   DG = [](nxn)<float> | The matrix of right-hand side coefficients
%
% Outputs:
%   delta_x = [](nx1)<float> | The answer with minimal residual (desc. above)
%   step_residual = [](1x1)<float> | The best answer's value of norm(g_i - G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Turn off warnings since some of these are often ill-conditioned
    warning('off', 'all');

    %% Solve the linear system using generalized min. residual method
    % Precondition DG to be better behaved, with 'P' row permutation,
    % 'R' row scaling, and 'C' column scaling.  Apply similar tx to 'G'
    [P, R, C] = equilibrate(DG);
    B = R * P * DG * C;
    d = R * P * G;

    % Try solving a conditioned problem
    [delta_y, ~, ~, ~, ~] = gmres(B, d, [], 1e-11, []);
    % Unconvert solved delta_y back into originally sought delta-x
    delta_x1 = C * delta_y;

    % Try solving with an unconditioned problem
    [delta_x2, ~, ~, ~, ~] = gmres(DG, G, [], 1e-11, []);

    %% Solve the linear system using Gauss-Jordan
    delta_x3 = rref([DG, G]);
    delta_x3 = delta_x3(:, end);

    %% Solve the linear system using mldivide (LU, Hessenberg, or Cholesky)
    % Warnings are suppressed to avoid annoying msg. on numerical
    % bad conditioning.  If condition number is poor, this method
    % won't be chosen over Gauss-Jordan... so I don't care if it
    % gives a shitty answer!
    delta_x4 = DG \ G;

    %% Determine which solution works the best
    delta_x_pos = [delta_x1, delta_x2, delta_x3, delta_x4];

    [step_residual, I] = min(vecnorm(DG*delta_x_pos-G, 2, 1));

    delta_x = delta_x_pos(:, I);

    % Re-enable warnings, mostly for the calling function
    warning('on', 'all');
end
