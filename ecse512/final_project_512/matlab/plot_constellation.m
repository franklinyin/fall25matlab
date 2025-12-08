function plot_constellation(z, S, ttl)
% Plot constellation scatter
    plot(real(z), imag(z), '.', 'MarkerSize', 6);
    hold on;
    grid on;
    axis equal;
    plot(real(S), imag(S), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('Re\{ \cdot \}');
    ylabel('Im\{ \cdot \}');
    title(ttl);
    legend('Samples','4-QAM ref','Location','best');
end
