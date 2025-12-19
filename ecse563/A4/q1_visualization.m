function q1_visualization(loads, c0, a, b, gmin, gmax, toler)
% helper function for production visualization in part 1

fprintf('Production vidualization\n');
figure('Name', 'q1 - production vidualization', 'Position', [100 100 1200 800]);

for idx = 1:length(loads)
    subplot(2, 2, idx);
    d = loads(idx);
    
    % Call ED to get results for this load
    [g, ~, lam] = ed(c0, a, b, gmin, gmax, d, toler);
    
    % Plot marginal cost curves for each generator
    g_range1 = linspace(gmin(1), gmax(1), 100);
    g_range2 = linspace(gmin(2), gmax(2), 100);
    g_range3 = linspace(gmin(3), gmax(3), 100);
    
    mc1 = a(1) + b(1) * g_range1;
    mc2 = a(2) + b(2) * g_range2;
    mc3 = a(3) + b(3) * g_range3;
    
    hold on;
    plot(g_range1, mc1, 'b-', 'LineWidth', 2, 'DisplayName', 'Gen 1 MC curve');
    plot(g_range2, mc2, 'r-', 'LineWidth', 2, 'DisplayName', 'Gen 2 MC curve');
    plot(g_range3, mc3, 'g-', 'LineWidth', 2, 'DisplayName', 'Gen 3 MC curve');
    
    % Plot lambda line
    yline(lam, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('\\lambda = %.3f', lam));
    
    % Plot actual dispatch points
    plot(g(1), a(1) + b(1)*g(1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
        'DisplayName', sprintf('Gen 1: %.2f MW', g(1)));
    plot(g(2), a(2) + b(2)*g(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('Gen 2: %.2f MW', g(2)));
    plot(g(3), a(3) + b(3)*g(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', ...
        'DisplayName', sprintf('Gen 3: %.2f MW', g(3)));
    
    grid on;
    xlabel('Generation Output (MW)', 'FontSize', 11);
    ylabel('Marginal Cost ($/MWh)', 'FontSize', 11);
    title(sprintf('d = %d MW, \\lambda = %.3f $/MWh', d, lam), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 9);
    hold off;
end
sgtitle('q1 - production vidualization', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'q1_production_vidualization.png');

end

