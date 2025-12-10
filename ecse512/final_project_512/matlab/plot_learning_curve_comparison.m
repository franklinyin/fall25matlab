function plot_learning_curve_comparison(cfg, channels, SNRdBs, resultsDir, comparisonType)
% Helper function to generate learning curve comparison plots
%
% Inputs:
%   cfg: base configuration struct
%   channels: cell array of channel types to compare (or single channel if comparing SNRs)
%   SNRdBs: array of SNR values to compare (or single SNR if comparing channels)
%   resultsDir: directory to save results
%   comparisonType: 'snr' or 'channel' to determine what we're comparing

    [S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));
    
    figure('Name', 'Learning Curve Comparison', 'Color', 'w', 'Position', [100 100 800 500]);
    hold on;
    
    legendEntries = {};
    colors = lines(max(length(channels), length(SNRdBs)));
    
    if strcmp(comparisonType, 'snr')
        % Comparing different SNRs for a single channel
        channelType = channels{1};
        cfg.channel.type = channelType;
        [h, hDescr] = generate_channel(cfg.channel);
        h = normalize_channel(h);
        D = pick_decision_delay(h, cfg.equalizerLenN);
        
        for idx = 1:length(SNRdBs)
            snr = SNRdBs(idx);
            [~, ~, ~, ~, ~, e, ~, ~, ~, mse, idxTrainEnd] = ...
                simulate_one_run(cfg, h, D, snr);
            
            plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, ...
                'Color', colors(idx,:));
            legendEntries{idx} = sprintf('SNR = %d dB', round(snr));
        end
        
        % Add training/decision switch line (same for all)
        xline(idxTrainEnd, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
        
        xlabel('Symbol index n', 'FontSize', 11);
        ylabel('MSE (dB)', 'FontSize', 11);
        title(sprintf('Learning Curves for %s Channel at Different SNRs', hDescr), 'FontSize', 12);
        legend(legendEntries, 'Location', 'best', 'FontSize', 10);
        grid on;
        
        filename = sprintf('learning_curve_%s_SNR_comparison.png', channelType);
        
    elseif strcmp(comparisonType, 'channel')
        % Comparing different channels at a single SNR
        snr = SNRdBs(1);
        
        for idx = 1:length(channels)
            channelType = channels{idx};
            cfg.channel.type = channelType;
            [h, hDescr] = generate_channel(cfg.channel);
            h = normalize_channel(h);
            D = pick_decision_delay(h, cfg.equalizerLenN);
            
            [~, ~, ~, ~, ~, e, ~, ~, ~, mse, idxTrainEnd] = ...
                simulate_one_run(cfg, h, D, snr);
            
            plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, ...
                'Color', colors(idx,:));
            legendEntries{idx} = sprintf('%s channel', hDescr);
        end
        
        % Add training/decision switch line (same for all)
        xline(idxTrainEnd, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
        
        xlabel('Symbol index n', 'FontSize', 11);
        ylabel('MSE (dB)', 'FontSize', 11);
        title(sprintf('Learning Curves for Different Channels at SNR = %d dB', round(snr)), 'FontSize', 12);
        legend(legendEntries, 'Location', 'best', 'FontSize', 10);
        grid on;
        
        filename = sprintf('learning_curve_channel_comparison_SNR%ddB.png', round(snr));
    end
    
    saveas(gcf, fullfile(resultsDir, filename));
    fprintf('Saved: %s\n', filename);
end
