function q1_min_output(loads, c0, a, b, gmin, gmax, toler)
% helper function for Fixed Cost Recovery Analysis for Economic Dispatch

fprintf('Fixed Cost Recovery Analysis:\n');

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

