function [x] = CVaR_robust(rets)
%CVAR_ROBUST Summary of this function goes here
%   Detailed explanation goes here
    r_alpha = 0.9;
    alpha = 0.95;

    options = optimoptions('fmincon', 'TolFun', 1e-9, 'MaxFunctionEvaluations', 100000);

    n = size(rets,2);
    N = size(rets, 1);
    S = N; 
    mu = (geomean(rets + 1) - 1)';
    Q = cov(rets);
    lambda = 0.1;

    % Defining bounds
    lb = [];
    ub = [];

    % Define the inequality constraint matrices A and b
    A = -[zeros(S, n) eye(S) zeros(S, 1);rets eye(S) ones(S,1)];
    b = [zeros(S, 1); zeros(S, 1)];

    % Define the equality constraint matrices A_eq and b_eq
    Aeq = [ones(1,n) zeros(1,S) 0];
    beq = 1;
    
    % Uncertainty set size
    Theta = diag(Q)./N;

    % Square root of Theta
    sqrtTh = sqrt(Theta); 

    % Scaling parameter epsilon for uncertainty set
    ep = sqrt(chi2inv(r_alpha,n));
    
    x0 = [(1/n)*ones(n, 1); ones(S, 1); 1];

    x = fmincon(@(x)objFun(x, mu, Q, lambda, sqrtTh, ep), x0, A, b, ...
                Aeq, beq, lb, ub, @(x)nonlcon(x), options);
    x = x(1:20);
    
end
function f = objFun(x, mu, Q, lambda, sqrtTh, ep)
  k = (1 / ((1 - 0.95) * 5000));
  f = x(end) + k*sum(x(21:end-1)) - lambda.*mu'*x(1:20)+ ep*norm(sqrtTh.*x(1:20));
end

function [c,ceq] = nonlcon(x)

    c = [];
    ceq = [];

end
