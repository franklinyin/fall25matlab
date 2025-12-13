function [g, C, lambda] = ed(co, a, b, gmin, gmax, d, toler)

    % Ensure column vectors
    co   = co(:);  a   = a(:);  b   = b(:);
    gmin = gmin(:); gmax = gmax(:);

    N = length(a);
    if any([length(co) length(b) length(gmin) length(gmax)] ~= N)
        error('All generator parameter vectors must have the same length.');
    end

    % Check feasibility of demand with limits
    if d < sum(gmin) || d > sum(gmax)
        error('Demand d is outside feasible range [sum(gmin), sum(gmax)].');
    end

    % Incremental costs at limits (for lambda bracketing)
    mc_min = a + b .* gmin;   % at gmin
    mc_max = a + b .* gmax;   % at gmax

    lambda_low  = min(mc_min);
    lambda_high = max(mc_max);

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
    C = sum(co + a .* g + 0.5 * b .* (g.^2));
end
