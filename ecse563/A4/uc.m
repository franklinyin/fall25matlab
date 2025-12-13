
function [u, g, C, lambda] = uc(c0, a, b, gmin, gmax, d, toler)
%UC Static unit commitment by full enumeration + ED on committed units
%   Inputs:
%     c0,a,b : column vectors (size N) of cost coefficients
%     gmin,gmax : column vectors (size N) of generator limits
%     d : total demand (scalar)
%     toler : power balance tolerance (scalar, e.g., 0.5 MW)
%   Outputs:
%     u : commitment status vector (0/1)
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

N = length(c0);
if any([length(a) length(b) length(gmin) length(gmax)] ~= N)
    error('All generator parameter vectors must have the same length.');
end

nComb = 2^N;

C_opt = inf;
u_opt = [];
g_opt = [];
lambda_opt = NaN;

% Enumerate all possible on/off combinations
for idx = 0:(nComb-1)
    % Binary representation -> commitment vector u in {0,1}^N
    bits = dec2bin(idx, N) - '0';
    u = bits(:);
    
    % At least one unit must be on if d>0
    if d > 0 && sum(u) == 0
        continue;
    end
    
    % Effective limits and fixed costs
    gmin_u = u .* gmin;
    gmax_u = u .* gmax;
    c0_u = u .* c0;
    
    % Solve ED for this commitment (units with u=0 forced at 0 MW)
    try
        [g, C_ed, lambda] = ed(c0_u, a, b, gmin_u, gmax_u, d, toler);
    catch
        % ED might throw if infeasible; skip this combination
        continue;
    end
    
    % Keep best (lowest-cost) feasible commitment
    if C_ed < C_opt
        C_opt = C_ed;
        u_opt = u;
        g_opt = g;
        lambda_opt = lambda;
    end
end

if isempty(u_opt)
    error('UC:NoFeasible', 'No feasible unit commitment found for this demand.');
end

% Return optimal solution
u = u_opt;
g = g_opt;
C = C_opt;
lambda = lambda_opt;
end
