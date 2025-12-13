function [u_opt, g_opt, C_opt, lambda_opt] = uc(co, a, b, gmin, gmax, d, toler)

    % Ensure column vectors
    co   = co(:);  a   = a(:);  b   = b(:);
    gmin = gmin(:); gmax = gmax(:);

    N = length(co);
    if any([length(a) length(b) length(gmin) length(gmax)] ~= N)
        error('All generator parameter vectors must have the same length.');
    end

    nComb = 2^N;

    C_opt      = inf;
    u_opt      = [];
    g_opt      = [];
    lambda_opt = NaN;

    % Enumerate all possible on/off combinations
    for idx = 0:(nComb-1)
        % Binary representation -> commitment vector u in {0,1}^N
        bits = dec2bin(idx, N) - '0' ;  % row, MSB first --> 000 to 111 for example each iteration
        u = bits(:);                    % make it column

        % At least one unit must be on if d>0
        if d > 0 && sum(u) == 0
            continue;
        end

        % Effective limits and fixed costs
        gmin_u = u .* gmin;
        gmax_u = u .* gmax;
        co_u   = u .* co;

        % Solve ED for this commitment (units with u=0 forced at 0 MW)
        try
            [g, C_ed, lambda] = ed(co_u, a, b, gmin_u, gmax_u, d, toler);
        catch
            % ED might throw if infeasible; skip this combination
            continue;
        end

        % Keep best (lowest-cost) feasible commitment
        if C_ed < C_opt
            C_opt      = C_ed;
            u_opt      = u;
            g_opt      = g;
            lambda_opt = lambda;
        end
    end

    if isempty(u_opt)
        error('UC:NoFeasible', 'No feasible unit commitment found for this demand.');
    end
end
