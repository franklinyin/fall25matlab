
function [u, g, C, lambda] = uc(c0, a, b, gmin, gmax, d, toler)
%UC Static unit commitment by full enumeration + ED on committed units
%   Inputs/outputs as in ED, plus u ∈ {0,1}^N (on/off statuses).
%   Feasible set requires sum(gmin(u)) <= d <= sum(gmax(u)).
%
%   Reference: ECSE 563 notes (Short‑Term Generation Optimization).

N = numel(a);
bestC = inf;
best = struct('u', [], 'g', [], 'C', inf, 'lambda', NaN);

for mask = 1:(2^N - 1)
    u_try = bitget(mask, 1:N)';  % binary vector
    gmin_u = gmin .* u_try;
    gmax_u = gmax .* u_try;
    if d < sum(gmin_u) - 1e-9 || d > sum(gmax_u) + 1e-9
        continue; % infeasible combination
    end
    % Run ED on committed set (off units are fixed at 0 via gmin=gmax=0)
    [g_try, C_try, lam_try] = ed(c0, a, b, gmin_u, gmax_u, d, toler);
    % Add fixed/no-load costs only for units that are on
    C_try = sum(c0 .* u_try) + sum(a.*g_try + 0.5*b.*(g_try.^2));
    if C_try < bestC - 1e-8
        bestC = C_try;
        best.u = u_try; best.g = g_try; best.C = C_try; best.lambda = lam_try;
    end
end

if isfinite(bestC)
    u = best.u; g = best.g; C = best.C; lambda = best.lambda;
else
    error('No feasible unit commitment.');
end
end
