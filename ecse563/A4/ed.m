
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
%
%   Reference: ECSE 563 notes (Short-Term Generation Optimization).

% Ensure column vectors
c0 = c0(:); 
a = a(:); 
b = b(:);
gmin = gmin(:); 
gmax = gmax(:);

N = length(a);
if any([length(c0) length(b) length(gmin) length(gmax)] ~= N)
    error('All generator parameter vectors must have the same length.');
end

% Feasibility check
if d < sum(gmin) || d > sum(gmax)
    error('Demand d is outside feasible range [sum(gmin), sum(gmax)].');
end

% Incremental costs at limits (for lambda bracketing)
mc_min = a + b .* gmin;
mc_max = a + b .* gmax;

lambda_low = min(mc_min);
lambda_high = max(mc_max);

% Lambda-iteration using bisection
maxiter = 100;
for k = 1:maxiter
    lambda = 0.5 * (lambda_low + lambda_high);
    
    % Unconstrained dispatch for this lambda
    g = (lambda - a) ./ b;
    
    % Enforce generator limits
    g = max(gmin, min(gmax, g));
    
    mismatch = sum(g) - d;
    
    if abs(mismatch) <= toler
        break;
    end
    
    if mismatch > 0
        % Too much generation -> lambda too high
        lambda_high = lambda;
    else
        % Not enough generation -> lambda too low
        lambda_low = lambda;
    end
end

if abs(sum(g) - d) > toler
    warning('ED:NoConverge', ...
        'Lambda-iteration reached maxiter without meeting tolerance.');
end

% Total cost
C = sum(c0 + a .* g + 0.5 * b .* (g.^2));
end
