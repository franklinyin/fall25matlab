function plot_ser_snr_comparison(cfg, channels, SNRdB_grid, resultsDir)
% Generate SER vs SNR comparison for multiple channels.

    figure;
    hold on;
    cols = lines(length(channels));
    legs = {};
    
    for ic = 1:length(channels)
        chtype = channels{ic};
        cfg.channel.type = chtype;
        [h, hdesc] = generate_channel(cfg.channel);
        h = normalize_channel(h);
        D = pick_decision_delay(h, cfg.equalizerLenN);
        
        SER = zeros(size(SNRdB_grid));
        for is = 1:numel(SNRdB_grid)
            sdb = SNRdB_grid(is);
            smc = zeros(1, cfg.numMC);
            for it = 1:cfg.numMC
                [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, sval] = simulate_one_run(cfg, h, D, sdb);
                smc(it) = sval;
            end
            SER(is) = mean(smc);
        end
        
        semilogy(SNRdB_grid, SER, '-o', 'LineWidth', 1.5, 'Color', cols(ic,:), ...
            'MarkerSize', 8, 'MarkerFaceColor', cols(ic,:));
        legs{ic} = sprintf('%s channel', hdesc);
    end
    
    grid on;
    xlabel('SNR (dB)', 'FontSize', 11);
    ylabel('Symbol Error Rate (SER)', 'FontSize', 11);
    title(sprintf('SER vs SNR for Different Channels (N=%d, \\mu=%.4f)', ...
        cfg.equalizerLenN, cfg.mu), 'FontSize', 12);
    legend(legs, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
    hold off;
    
    fname = fullfile(resultsDir, 'SER_vs_SNR_channel_comparison.png');
    print(gcf, fname, '-dpng', '-r300');
end







