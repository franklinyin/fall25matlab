% Q1 implementation
function [g, C, lambda] = ed(c0, a, b, gmin, gmax, d, toler)
%  economic dispatch by lambda-iteration (quadratic costs)
    %   Inputs:
    %     c0,a,b : column vectors (size N) of cost coefficients
    %     gmin,gmax : column vectors (size N) of generator limits
    %     d : total demand (scalar)
    %     toler : power balance tolerance (scalar, e.g., 0.5 MW)
    %   Outputs:
    %     g : dispatch vector (MW)
    %     C : total cost ($/h)
    %     lambda : power balance Lagrange multiplier ($/MWh)


    N = length(a);
    
    % Feasibility check
    if d < sum(gmin) || d > sum(gmax)
        error('Demand d is outside feasible range [sum(gmin), sum(gmax)].');
    end
    
    % -------- lambda-iteration algorithm (gradient method) --------
    
    % Incremental costs at limits (C'_i(g))
    mc_min = a + b .* gmin;
    mc_max = a + b .* gmax;
    
    % Initial lambda (from notes):
    % lambda^0 = ( d + sum (a_i / b_i) ) / sum (1 / b_i)
    lambda = (d + sum(a ./ b)) / sum(1 ./ b);
    
    % Parameters
    beta    = 0.03;      % step size
    maxiter = 10000;     % safety cap on iterations
    
    % Initialization
    Delta  = Inf;        % power mismatch
    g      = zeros(N,1); % dispatch vector
    
    k = 0;
    while abs(Delta) > toler && k < maxiter
        k = k + 1;
        lambda_old = lambda;   % lambda used to compute this iteration's dispatch
        
        % Compute g_i^k from lambda^k with limit checks
        for i = 1:N
            if mc_min(i) >= lambda_old
                g(i) = gmin(i);
            elseif mc_max(i) <= lambda_old
                g(i) = gmax(i);
            else
                g(i) = (lambda_old - a(i)) / b(i);
            end
        end
        
        % Power balance mismatch Delta^k = sum g_i^k - d
        Delta = sum(g) - d;

        lambda = lambda_old - beta * Delta;
    end
    
    % Use the lambda that produced the final dispatch g
    lambda = lambda_old;
    
    if abs(Delta) > toler % for part uc on feasibility
        warning('No Convergence');
    end
    
    C = sum(c0 + a .* g + 0.5 * b .* (g.^2));
end
