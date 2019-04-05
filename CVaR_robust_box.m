function [x] = CVaR_robust_box(rets)
%CVAR_ROBUST Summary of this function goes here
%   Detailed explanation goes here
    r_alpha = 0.9;
    alpha = 0.95;

    options = optimoptions('fmincon', 'TolFun', 1e-9, 'MaxFunctionEvaluations', 30000);

    n = size(rets, 2);
    N = size(rets, 1);
    S = N; 
    mu = (geomean(rets + 1) - 1)';
    Q = cov(rets);
    lambda = 0.1;

    % Defining bounds
    lb = [];
    ub = [];

    % Define the inequality constraint matrices A and b
    A = -[zeros(S, n) zeros(S, n) eye(S) zeros(S, 1); 
    rets zeros(S, n) eye(S) ones(S, 1);
    eye(n) eye(n) zeros(n, S) zeros(n, 1);
    -eye(n) eye(n) zeros(n, S) zeros(n, 1)];
    
    b = [zeros(S, 1); zeros(S, 1); zeros(n, 1); zeros(n, 1)];

    % Define the equality constraint matrices A_eq and b_eq
    Aeq = [ones(1, n) zeros(1, n) zeros(1, S) 0];
    beq = 1;
    
    % Uncertainty set size
    Theta = diag(Q) ./ N;

    % Square root of Theta
    sqrtTh = sqrt(Theta); 

    % Scaling parameter epsilon for uncertainty set
    ep = norminv(1 - (1 - r_alpha) / 2, 0, 1);
    
    k = (1 / ((1 - alpha) * S));

    
    c = [-lamda * mu;
         lambda * ep .* sqrtTh' * ones(S,1);
          k * ones(S,1);
          1       ];
    y = linprog( c, A, b, Aeq, beq, lb, ub );
 

    x = y(1:n);
    
end

