function [x] = CVaR_robust_box(rets)
%CVAR_ROBUST Summary of this function goes here
%   Detailed explanation goes here
% rets = [0.0199857690907118,-0.0256079854742501;
%     0.0411098790489003,0.0167472169476907;
%     0.0411098790489003,0.0167472169476907]
    r_alpha = 0.9;
    alpha = 0.95;
    
    n = size(rets, 2);
    N = size(rets, 1);
    S = N; 
    mu = (geomean(rets + 1) - 1)';
    Q = cov(rets);
    lambda = 0.1;

    % Defining bounds
    lb = [-inf(n, 1); -inf(n, 1); zeros(S, 1); -inf];
    ub = [];

    % Define the inequality constraint matrices A and b
    A = [zeros(1, n) zeros(1, n) -ones(1, S) zeros(1, 1);
        -rets zeros(size(rets)) -eye(S) -ones(S, 1);
        -eye(n) -eye(n) zeros(n, S) zeros(n, 1);
        eye(n)  -eye(n) zeros(n, S) zeros(n, 1)];
    b = zeros(S+ 1 + n + n, 1);
    % Define the equality constraint matrices A_eq and b_eq
    Aeq = [ones(1, n) zeros(1, n) zeros(1, S) 0];
    beq = 1;
    
    % Uncertainty set size
    Theta = diag(Q) ./ N;

    % Square root of Theta
    sqrtTh = sqrt(Theta); 

    % Scaling parameter epsilon for uncertainty set
    ep = norminv(1 - (1 - r_alpha) / 2, 0, 1);
    
    k = (1 / ((1 - alpha) * S))


    c = [-lambda .* mu;
         lambda * ep .* sqrtTh;
          k .* ones(S,1);
          1       ];
    
    x = linprog( c, A, b, Aeq, beq, lb, ub );
 
    x = x(1:n);
    
end

