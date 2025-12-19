% Q2
function [u, g, C, lambda] = uc(c0, a, b, gmin, gmax, d, toler)
    % uc - unit commitment by enumeration + ED
    
    N = numel(c0);
    % if any([numel(a) numel(b) numel(gmin) numel(gmax)] ~= N)
    %     error('All generator parameter vectors must have the same length.');
    % end
    
    ncomb = 2^N;
    
    Cbest = inf;
    ubest = [];
    gbest = [];
    lambest = NaN;
    
    for idx = 0:(ncomb-1)
        bits = dec2bin(idx, N) - '0';
        u = bits(:);
        
        if d > 0 && sum(u) == 0
            continue;
        end
        
        gmin_u = u .* gmin;
        gmax_u = u .* gmax;
        c0_u = u .* c0;
        
        try
            [g, C_ed, lambda] = ed(c0_u, a, b, gmin_u, gmax_u, d, toler);
        catch
            continue;
        end
        
        if C_ed < Cbest
            Cbest = C_ed;
            ubest = u;
            gbest = g;
            lambest = lambda;
        end
    end
    
    if isempty(ubest)
        error('No UC found');
    end
    
    u = ubest;
    g = gbest;
    C = Cbest;
    lambda = lambest;
end
