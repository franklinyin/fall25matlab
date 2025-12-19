function q2_q1_profit_comparison(loads, c0, a, b, gmin, gmax, toler)
% Q2_Q1_PROFIT_COMPARISON Profit Comparison: Unit Commitment vs Economic Dispatch
% 
% Inputs:
%   loads  - Array of load values (MW)
%   c0     - Fixed cost coefficients (vector)
%   a      - Linear cost coefficients (vector)
%   b      - Quadratic cost coefficients (vector)
%   gmin   - Minimum generation limits (vector)
%   gmax   - Maximum generation limits (vector)
%   toler  - Tolerance for ED/UC convergence
%
% Compares profits between UC and ED methods when generators are remunerated
% at marginal cost lambda.

fprintf('\n----------------------------------------------------------------------\n');
fprintf('PROFIT COMPARISON: Unit Commitment vs Economic Dispatch\n');
fprintf('----------------------------------------------------------------------\n');
fprintf('Generators remunerated at marginal cost lambda\n');
fprintf('Profit_i = lambda * g_i - Cost_i(g_i)\n\n');

fprintf('%-8s', 'Load');
fprintf('%-12s', 'Method');
fprintf('%-18s', 'Total Profit');
fprintf('%-18s', 'Commitment');
fprintf('%-18s\n', 'Cost');
fprintf('%-8s', '(MW)');
fprintf('%-12s', '');
fprintf('%-18s', '($/h)');
fprintf('%-18s', '');
fprintf('%-18s\n', '($/h)');
fprintf('----------------------------------------------------------------------\n');

for idx = 1:length(loads)
    d = loads(idx);
    
    % Call ED to get results
    [g_ed, C_ed, lam_ed] = ed(c0, a, b, gmin, gmax, d, toler);
    
    % Calculate ED profits (all units committed: u = [1,1,1])
    prof_ed = lam_ed * g_ed - (c0 + a .* g_ed + 0.5 * b .* (g_ed.^2));
    total_profit_ed = sum(prof_ed);
    
    % Call UC to get results
    [u_uc, g_uc, C_uc, lam_uc] = uc(c0, a, b, gmin, gmax, d, toler);
    
    % Calculate UC profits
    prof_uc = lam_uc * g_uc - (c0.*u_uc + a.*g_uc + 0.5*b.*g_uc.^2);
    total_profit_uc = sum(prof_uc);
    
    % ED row
    fprintf('%-8.0f %-12s %-18.2f %-18s %-18.2f\n', ...
        d, 'ED', total_profit_ed, '[1,1,1]', C_ed);
    
    % UC row
    u_str = sprintf('[%d,%d,%d]', u_uc(1), u_uc(2), u_uc(3));
    fprintf('%-8s %-12s %-18.2f %-18s %-18.2f\n', ...
        '', 'UC', total_profit_uc, u_str, C_uc);
    
    % Difference
    profit_diff = total_profit_uc - total_profit_ed;
    fprintf('%-8s %-12s %-18.2f\n\n', '', 'UC - ED', profit_diff);
end

end

