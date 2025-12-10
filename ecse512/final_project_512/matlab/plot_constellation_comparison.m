function plot_constellation_comparison(cfg, h, hDescr, D, SNRdB, resultsDir)
% Generate and save constellation diagram for given SNR and channel.

    [S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));
    [~, y, ~, ~, xhat, ~, ~, ~, ~, ~, ~] = simulate_one_run(cfg, h, D, SNRdB);
    
    % extract channel type for filename
    chlabel = lower(strrep(hDescr, '‑', '-'));
    if contains(chlabel, 'low')
        chtype = 'low';
        chname = 'Low ISI';
    elseif contains(chlabel, 'severe')
        chtype = 'severe';
        chname = 'Severe ISI';
    elseif contains(chlabel, 'mild')
        chtype = 'mild';
        chname = 'Mild ISI';
    else
        chtype = lower(strrep(chlabel, ' ', '_'));
        chname = hDescr;
    end
    
    figure('Name', sprintf('Constellations %s SNR%ddB', hDescr, round(SNRdB)), 'Color', 'w');
    
    % left: input
    subplot(1,2,1);
    plot(real(y), imag(y), '.', 'MarkerSize', 6); hold on; grid on; axis equal;
    plot(real(S), imag(S), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('Re\{ \cdot \}'); ylabel('Im\{ \cdot \}'); title('Rx input y[n]');
    legend('Samples','4‑QAM ref','Location','best');
    
    % right: output
    subplot(1,2,2);
    plot(real(xhat), imag(xhat), '.', 'MarkerSize', 6); hold on; grid on; axis equal;
    plot(real(S), imag(S), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('Re\{ \cdot \}'); ylabel('Im\{ \cdot \}'); title('Equalizer output \hat{x}[n]');
    legend('Samples','4‑QAM ref','Location','best');
    
    sgtitle(sprintf('%s Channel at SNR = %d dB', chname, round(SNRdB)), 'FontSize', 14, 'FontWeight', 'bold');
    
    fname = sprintf('constellation_%s_SNR%ddB.png', chtype, round(SNRdB));
    saveas(gcf, fullfile(resultsDir, fname));
end

