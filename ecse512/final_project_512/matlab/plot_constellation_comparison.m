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
    subplot(1,2,1);
    plot_constellation(y, S, 'Rx input y[n]');
    subplot(1,2,2);
    plot_constellation(xhat, S, 'Equalizer output \hat{x}[n]');
    sgtitle(sprintf('%s Channel at SNR = %d dB', chname, round(SNRdB)), 'FontSize', 14, 'FontWeight', 'bold');
    
    fname = sprintf('constellation_%s_SNR%ddB.png', chtype, round(SNRdB));
    saveas(gcf, fullfile(resultsDir, fname));
    fprintf('Saved: %s\n', fname);
end
