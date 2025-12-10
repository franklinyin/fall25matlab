function plot_constellation_comparison(cfg, h, hDescr, D, SNRdB, resultsDir)
% Helper function to generate and save constellation diagram for a given configuration
%
% Inputs:
%   cfg: configuration struct
%   h: channel impulse response
%   hDescr: channel description string
%   D: decision delay
%   SNRdB: SNR in dB for this run
%   resultsDir: directory to save results

    % Get constellation reference
    [S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));
    
    % Run simulation
    [~, y, ~, ~, xhat, ~, ~, ~, ~, ~, ~] = simulate_one_run(cfg, h, D, SNRdB);
    
    % Extract just the channel type (low, mild, severe)
    channelLabel = lower(strrep(hDescr, '‑', '-'));
    if contains(channelLabel, 'low')
        channelType = 'low';
        channelName = 'Low ISI';
    elseif contains(channelLabel, 'severe')
        channelType = 'severe';
        channelName = 'Severe ISI';
    elseif contains(channelLabel, 'mild')
        channelType = 'mild';
        channelName = 'Mild ISI';
    else
        channelType = lower(strrep(channelLabel, ' ', '_'));
        channelName = hDescr;
    end
    
    % Create figure with big title
    figure('Name', sprintf('Constellations %s SNR%ddB', hDescr, round(SNRdB)), 'Color', 'w');
    
    subplot(1,2,1);
    plot_constellation(y, S, 'Rx input y[n]');
    subplot(1,2,2);
    plot_constellation(xhat, S, 'Equalizer output \hat{x}[n]');
    
    % Add overall title
    sgtitle(sprintf('%s Channel at SNR = %d dB', channelName, round(SNRdB)), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save with descriptive filename
    filename = sprintf('constellation_%s_SNR%ddB.png', channelType, round(SNRdB));
    saveas(gcf, fullfile(resultsDir, filename));
    fprintf('Saved: %s\n', filename);
end
