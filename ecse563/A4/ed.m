% Q1 implementation
function [g, C, lambda] = ed(c0, a, b, gmin, gmax, d, toler)
%ED Economic dispatch by lambda-iteration (quadratic costs)
%   Cost Ci(gi) = c0_i + a_i*gi + 0.5*b_i*gi^2
%   Inputs:
%     c0,a,b : column vectors (size N) of cost coefficients
%     gmin,gmax : column vectors (size N) of generator limits
%     d : total demand (scalar)
%     toler : power balance tolerance (scalar, e.g., 0.5 MW)
%   Outputs:
%     g : dispatch vector (MW)
%     C : total cost ($/h)
%     lambda : power balance Lagrange multiplier ($/MWh)

% Ensure column vectors
% c0   = c0(:);
% a    = a(:);
% b    = b(:);
% gmin = gmin(:);
% gmax = gmax(:);

N = length(a);
% if any([length(c0) length(b) length(gmin) length(gmax)] ~= N)
%     error('All generator parameter vectors must have the same length.');
% end

% Feasibility check
if d < sum(gmin) || d > sum(gmax)
    error('Demand d is outside feasible range [sum(gmin), sum(gmax)].');
end

% -------- λ-iteration algorithm (gradient method) --------

% Incremental costs at limits (C'_i(g))
mc_min = a + b .* gmin;
mc_max = a + b .* gmax;

% Initial lambda (from notes):
% λ^0 = ( d + Σ (a_i / b_i) ) / Σ (1 / b_i)
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
    lambda_old = lambda;   % λ used to compute this iteration's dispatch
    
    % Compute g_i^k from λ^k with limit checks
    for i = 1:N
        if mc_min(i) >= lambda_old
            g(i) = gmin(i);
        elseif mc_max(i) <= lambda_old
            g(i) = gmax(i);
        else
            g(i) = (lambda_old - a(i)) / b(i);
        end
    end
    
    % Power balance mismatch Δ^k = Σ g_i^k - d
    Delta = sum(g) - d;
    
    % Update λ for next iteration: λ^{k+1} = λ^k - β Δ^k
    lambda = lambda_old - beta * Delta;
end

% Use the λ that produced the final dispatch g
lambda = lambda_old;

if abs(Delta) > toler
    warning('ED:NoConverge', ...
        'Lambda-iteration reached maxiter without meeting tolerance.');
end

% Total cost
C = sum(c0 + a .* g + 0.5 * b .* (g.^2));
end
