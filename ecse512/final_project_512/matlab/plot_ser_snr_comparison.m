function plot_ser_snr_comparison(cfg, channels, SNRdB_grid, resultsDir)
% Helper function to generate SER vs SNR comparison plot for multiple channels
%
% Inputs:
%   cfg: base configuration struct
%   channels: cell array of channel types to compare
%   SNRdB_grid: array of SNR values to test
%   resultsDir: directory to save results

    figure('Name', 'SER vs SNR Comparison', 'Color', 'w', 'Position', [100 100 800 500]);
    hold on;
    
    colors = lines(length(channels));
    legendEntries = {};
    
    for chIdx = 1:length(channels)
        channelType = channels{chIdx};
        cfg.channel.type = channelType;
        
        % Generate and normalize channel
        [h, hDescr] = generate_channel(cfg.channel);
        h = normalize_channel(h);
        D = pick_decision_delay(h, cfg.equalizerLenN);
        
        fprintf('\nComputing SER vs SNR for %s channel...\n', hDescr);
        
        % Compute SER for each SNR
        SER = zeros(size(SNRdB_grid));
        for iS = 1:numel(SNRdB_grid)
            snrdb = SNRdB_grid(iS);
            ser_mc = zeros(1, cfg.numMC);
            
            for it = 1:cfg.numMC
                [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
                    simulate_one_run(cfg, h, D, snrdb);
                ser_mc(it) = ser_val;
            end
            
            SER(iS) = mean(ser_mc);
            fprintf('  SNR = %2d dB --> SER = %.3e\n', snrdb, SER(iS));
        end
        
        % Plot
        semilogy(SNRdB_grid, SER, '-o', 'LineWidth', 1.5, 'Color', colors(chIdx,:), ...
            'MarkerSize', 8, 'MarkerFaceColor', colors(chIdx,:));
        
        legendEntries{chIdx} = sprintf('%s channel', hDescr);
    end
    
    grid on;
    xlabel('SNR (dB)', 'FontSize', 11);
    ylabel('Symbol Error Rate (SER)', 'FontSize', 11);
    title(sprintf('SER vs SNR for Different Channels (N=%d, \\mu=%.4f)', ...
        cfg.equalizerLenN, cfg.mu), 'FontSize', 12);
    
    % Create and configure legend - force it to be visible
    legend(legendEntries, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
    
    hold off;
    
    % Save using print for better legend rendering
    filename = fullfile(resultsDir, 'SER_vs_SNR_channel_comparison.png');
    print(gcf, filename, '-dpng', '-r300');
    fprintf('\nSaved: SER_vs_SNR_channel_comparison.png\n');
end
