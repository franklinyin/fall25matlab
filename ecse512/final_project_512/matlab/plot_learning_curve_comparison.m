function plot_learning_curve_comparison(cfg, channels, SNRdBs, resultsDir, comparisonType)
% Generate learning curve comparison (either across SNR or across channels).

    figure;
    hold on;
    
    legs = {};
    cols = lines(max(length(channels), length(SNRdBs)));
    
    if strcmp(comparisonType, 'snr')
        chtype = channels{1};
        cfg.channel.type = chtype;
        [h, hdesc] = generate_channel(cfg.channel);
        h = normalize_channel(h);
        D = pick_decision_delay(h, cfg.equalizerLenN);
        
        for i = 1:length(SNRdBs)
            snr = SNRdBs(i);
            [~, ~, ~, ~, ~, ~, ~, ~, ~, mse, itr] = simulate_one_run(cfg, h, D, snr);
            plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, 'Color', cols(i,:));
            legs{i} = sprintf('SNR = %d dB', round(snr));
        end
        xline(itr, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
        xlabel('Symbol index n', 'FontSize', 11);
        ylabel('MSE (dB)', 'FontSize', 11);
        title(sprintf('Learning Curves for %s Channel at Different SNRs', hdesc), 'FontSize', 12);
        legend(legs, 'Location', 'best', 'FontSize', 10);
        grid on;
        fname = sprintf('learning_curve_%s_SNR_comparison.png', chtype);
        
    elseif strcmp(comparisonType, 'channel')
        snr = SNRdBs(1);
        for i = 1:length(channels)
            chtype = channels{i};
            cfg.channel.type = chtype;
            [h, hdesc] = generate_channel(cfg.channel);
            h = normalize_channel(h);
            D = pick_decision_delay(h, cfg.equalizerLenN);
            [~, ~, ~, ~, ~, ~, ~, ~, ~, mse, itr] = simulate_one_run(cfg, h, D, snr);
            plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, 'Color', cols(i,:));
            legs{i} = sprintf('%s channel', hdesc);
        end
        xline(itr, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
        xlabel('Symbol index n', 'FontSize', 11);
        ylabel('MSE (dB)', 'FontSize', 11);
        title(sprintf('Learning Curves for Different Channels at SNR = %d dB', round(snr)), 'FontSize', 12);
        legend(legs, 'Location', 'best', 'FontSize', 10);
        grid on;
        fname = sprintf('learning_curve_channel_comparison_SNR%ddB.png', round(snr));
    end
    
    saveas(gcf, fullfile(resultsDir, fname));
end

