function q1_min_output(loads, c0, a, b, gmin, gmax, toler)
% Q1_MIN_OUTPUT Fixed Cost Recovery Analysis for Economic Dispatch
% 
% Inputs:
%   loads  - Array of load values (MW)
%   c0     - Fixed cost coefficients (vector)
%   a      - Linear cost coefficients (vector)
%   b      - Quadratic cost coefficients (vector)
%   gmin   - Minimum generation limits (vector)
%   gmax   - Maximum generation limits (vector)
%   toler  - Tolerance for ED convergence
%
% Analyzes whether generators can recover fixed costs at the dispatch lambda.

fprintf('c) Fixed Cost Recovery Analysis:\n');
fprintf('-------------------------------------------------\n');
fprintf('For generator i to recover fixed costs when remunerated at lambda*g_i:\n');
fprintf('  Require: lambda * g_i >= c0_i  =>  g_i >= c0_i / lambda\n\n');

fprintf('%-8s', 'Load');
fprintf('%-12s', 'Lambda');
fprintf('%-12s', 'Gen 1 Min');
fprintf('%-12s', 'Gen 2 Min');
fprintf('%-12s\n', 'Gen 3 Min');
fprintf('%-8s', '(MW)');
fprintf('%-12s', '($/MWh)');
fprintf('%-12s', '(MW)');
fprintf('%-12s', '(MW)');
fprintf('%-12s\n', '(MW)');
fprintf('------------------------------------------------\n');

for idx = 1:length(loads)
    d = loads(idx);
    
    % Call ED to get results for this load
    [g, ~, lam] = ed(c0, a, b, gmin, gmax, d, toler);
    
    g_min_recovery = c0 ./ lam;
    
    fprintf('%-8.0f', d);
    fprintf('%-12.3f', lam);
    fprintf('%-12.2f', g_min_recovery(1));
    fprintf('%-12.2f', g_min_recovery(2));
    fprintf('%-12.2f\n', g_min_recovery(3));
    
    % Check if actual dispatch meets recovery requirement
    fprintf('          Actual dispatch: [%.2f, %.2f, %.2f] MW\n', g(1), g(2), g(3));
    meets = g >= g_min_recovery;
    fprintf('          Meets recovery: [%s, %s, %s]\n\n', ...
        bool2str(meets(1)), bool2str(meets(2)), bool2str(meets(3)));
end

end

function str = bool2str(val)
    if val
        str = 'YES';
    else
        str = 'NO';
    end
end

