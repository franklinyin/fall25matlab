% Q2
function [u, g, C, lambda] = uc(c0, a, b, gmin, gmax, d, toler)
    % uc - unit commitment by enumeration + ED
    % Inputs:
    %   c0,a,b : column vectors (size N) of cost coefficients
    %   gmin,gmax : column vectors (size N) of generator limits
    %   d : total demand (scalar)
    %   toler : power balance tolerance
    % Outputs:
    %   u : commitment vector (0/1)
    %   g : dispatch vector (MW)
    %   C : total cost ($/h)
    %   lambda : power balance Lagrange multiplier ($/MWh)
    
    N = numel(c0);
    
    ncomb = 2^N;
    
    Cbest = inf;
    ubest = [];
    gbest = [];
    lambest = NaN;
    
    for idx = 0:(ncomb-1) % enumerate all possibilities
        bits = dec2bin(idx, N) - '0';
        u = bits(:);
        
        if d > 0 && sum(u) == 0 % at least one unit producing
            continue;
        end
        
        gmin_u = u .* gmin;
        gmax_u = u .* gmax;
        c0_u = u .* c0;
        
        try %trying ed and reject for infeasibility
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
