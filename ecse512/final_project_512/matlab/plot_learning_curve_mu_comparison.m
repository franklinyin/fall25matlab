function plot_learning_curve_mu_comparison(cfg, channel, SNRdB, muGrid, resultsDir)
% Compare learning curves for different step sizes mu.

    cfg.channel.type = channel;
    [h, hdesc] = generate_channel(cfg.channel);
    h = normalize_channel(h);
    D = pick_decision_delay(h, cfg.equalizerLenN);
    
    figure;
    hold on;
    cols = lines(length(muGrid));
    legs = {};
    
    for i = 1:length(muGrid)
        cfg2 = cfg;
        cfg2.mu = muGrid(i);
        [~, ~, ~, ~, ~, ~, ~, ~, ~, mse, itr] = simulate_one_run(cfg2, h, D, SNRdB);
        plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth', 1.5, 'Color', cols(i,:));
        legs{i} = sprintf('\\mu = %.4f', muGrid(i));
    end
    
    xline(itr, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
    xlim([0 4000]);
    xlabel('Symbol index n', 'FontSize', 11);
    ylabel('MSE (dB)', 'FontSize', 11);
    title(sprintf('Learning Curves for Different Step Sizes (SNR=%d dB, %s)', round(SNRdB), hdesc), 'FontSize', 12);
    legend(legs, 'Location', 'best', 'FontSize', 10);
    grid on;
    hold off;
    
    fname = sprintf('learning_curve_mu_comparison_%s_SNR%ddB.png', channel, round(SNRdB));
    saveas(gcf, fullfile(resultsDir, fname));
end







