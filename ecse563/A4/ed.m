
function [g, C, lambda] = ed(c0, a, b, gmin, gmax, d, toler)
%ED Economic dispatch by lambda-iteration (quadratic costs)
%   Cost Ci(gi) = c0_i + a_i*gi + 0.5*b_i*gi^2
%   Inputs:
%     c0,a,b : column vectors (size N) of cost coefficients
%     gmin,gmax : column vectors (size N) of generator limits
%     d : total demand (scalar)
%     toler : power balance tolerance (scalar, e.g., 0.5 MW)
%   Outputs:
%     g : dispatch vector
%     C : total cost ($/h)
%     lambda : power balance Lagrange multiplier ($/MWh)
%
%   Reference: ECSE 563 notes (Short‑Term Generation Optimization).

N = numel(a);
c0 = c0(:); a = a(:); b = b(:); gmin = gmin(:); gmax = gmax(:);

% Feasibility check
if d < sum(gmin) - 1e-9 || d > sum(gmax) + 1e-9
    error('Infeasible demand: outside sum(gmin)/sum(gmax).');
end

% Bracket lambda using marginal costs at the bounds
lam_lo = min(a + b.*gmin) - 1000;
lam_hi = max(a + b.*gmax) + 1000;

% Monotone bisection on sum(g(lambda)) - d
for it = 1:200
    lam = 0.5*(lam_lo + lam_hi);
    g = (lam - a)./b;
    g = min(max(g, gmin), gmax);
    s = sum(g);
    if s < d
        lam_lo = lam;
    else
        lam_hi = lam;
    end
    if abs(lam_hi - lam_lo) < 1e-9 || abs(s - d) <= toler
        break;
    end
end
lambda = 0.5*(lam_lo + lam_hi);
g = (lambda - a)./b;
g = min(max(g, gmin), gmax);

% Total cost
C = sum(c0 + a.*g + 0.5*b.*(g.^2));
end
