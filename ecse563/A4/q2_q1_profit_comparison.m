function q2_q1_profit_comparison(loads, c0, a, b, gmin, gmax, toler)
% Helper function for profit comparison between UC and ED
fprintf('\n----------------------------------------------------------------------\n');
fprintf('PROFIT COMPARISON: Unit Commitment vs Economic Dispatch\n');
fprintf('----------------------------------------------------------------------\n');

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
    
    % call ED to get results
    [g_ed, C_ed, lam_ed] = ed(c0, a, b, gmin, gmax, d, toler);
    
    % calculate ED profits (all units committed: u = [1,1,1])
    prof_ed = lam_ed * g_ed - (c0 + a .* g_ed + 0.5 * b .* (g_ed.^2));
    total_profit_ed = sum(prof_ed);
    
    % call UC to get results
    [u_uc, g_uc, C_uc, lam_uc] = uc(c0, a, b, gmin, gmax, d, toler);
    
    % calculate UC profits
    prof_uc = lam_uc * g_uc - (c0.*u_uc + a.*g_uc + 0.5*b.*g_uc.^2);
    total_profit_uc = sum(prof_uc);
    
    % ED row
    fprintf('%-8.0f %-12s %-18.2f %-18s %-18.2f\n', ...
        d, 'ED', total_profit_ed, '[1,1,1]', C_ed);
    
    % UC row
    u_str = sprintf('[%d,%d,%d]', u_uc(1), u_uc(2), u_uc(3));
    fprintf('%-8s %-12s %-18.2f %-18s %-18.2f\n', ...
        '', 'UC', total_profit_uc, u_str, C_uc);
    
    % difference
    profit_diff = total_profit_uc - total_profit_ed;
    fprintf('%-8s %-12s %-18.2f\n\n', '', 'UC - ED', profit_diff);
end

end

