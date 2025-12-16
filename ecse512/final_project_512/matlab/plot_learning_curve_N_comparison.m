function plot_learning_curve_N_comparison(cfg, channel, SNRdB, Ngrid, resultsDir)
% Compare learning curves for different equalizer lengths N.

    cfg.channel.type = channel;
    [h, hdesc] = generate_channel(cfg.channel);
    h = normalize_channel(h);
    
    figure;
    hold on;
    cols = lines(length(Ngrid));
    legs = {};
    
    for i = 1:length(Ngrid)
        cfg2 = cfg;
        cfg2.equalizerLenN = Ngrid(i);
        D = pick_decision_delay(h, cfg2.equalizerLenN);
        [~, ~, ~, ~, ~, ~, ~, ~, ~, mse, itr] = simulate_one_run(cfg2, h, D, SNRdB);
        plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, 'Color', cols(i,:));
        legs{i} = sprintf('N = %d', Ngrid(i));
    end
    
    xline(itr, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
    xlim([0 4000]);
    xlabel('Symbol index n', 'FontSize', 11);
    ylabel('MSE (dB)', 'FontSize', 11);
    title(sprintf('Learning Curves for Different Equalizer Lengths (SNR=%d dB, %s)', round(SNRdB), hdesc), 'FontSize', 12);
    legend(legs, 'Location', 'best', 'FontSize', 10);
    grid on;
    hold off;
    
    fname = sprintf('learning_curve_N_comparison_%s_SNR%ddB.png', channel, round(SNRdB));
    saveas(gcf, fullfile(resultsDir, fname));
end



