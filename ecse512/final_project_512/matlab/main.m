% main script

clear; clc; close all;

cfg = struct();

% simulation parameters
cfg.M = 4; % QAM order (fixed to 4‑QAM)
cfg.sigma_s2 = 1; % symbol power (E{|x|^2})
cfg.trainLen = 1500; % training symbols
cfg.dataLen = 15000; % data symbols (decision‑directed phase)
cfg.equalizerLenN = 9; % LMS equalizer length (N taps)
cfg.mu = 0.015; % LMS step size (tune per scenario)
cfg.numMC = 5; % Monte‑Carlo trials (increase for stats)
cfg.SNRdB_grid = 0:3:18; % SNR sweep for SER curves
cfg.saveFigs = true; % save figures to ../results
cfg.randomSeed = 42; % reproducibility
cfg.channel.type = 'severe'; % 'low' | 'mild' | 'severe'
cfg.plotOneRunSNRdB = 30; % SNR for which to show a detailed run
cfg.centerTap = []; % optional target decision delay D (0..N-1)

rng(cfg.randomSeed);

% make output dir
resultsDir = fullfile('..','results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

% reference 4‑QAM constellation (unit power scaled to cfg.sigma_s2)
[S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));  % S is 4 symbols

% channel definition and decision delay
[h, hDescr] = generate_channel(cfg.channel);
h = normalize_channel(h);  % unit energy
K = numel(h);
fprintf('Channel (%s): h = %s\n', hDescr, mat2str(h,4));

% choose decision delay D to align the equalizer target
if isempty(cfg.centerTap)
    D = pick_decision_delay(h, cfg.equalizerLenN);
else
    D = min(max(cfg.centerTap,0), cfg.equalizerLenN-1);
end
fprintf('Decision delay D = %d (equalizer length N = %d)\n', D, cfg.equalizerLenN);

% one illustrative run for learning curve
SNR_show = cfg.plotOneRunSNRdB;
fprintf('\n=== Single run at SNR = %.1f dB for learning curve ===\n', SNR_show);
[x, y, w, noiseVar, xhat, e, d_used, c_hist, dec, mse, idxTrainEnd] = ...
    simulate_one_run(cfg, h, D, SNR_show);

% learning curve
figure('Name','Learning curve (one run)','Color','w');
plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth',1.2);
hold on; yline(10*log10(mean(mse(end-100:end))), '--');
xline(idxTrainEnd,'--');
xlabel('n'); ylabel('E[|e[n]|^2] (dB)'); grid on;
title(sprintf('Learning curve, SNR=%.1f dB, N=%d, \\mu=%.4f, D=%d', ...
      SNR_show, cfg.equalizerLenN, cfg.mu, D));
legend('MSE (smoothed)','Steady‑state (avg last 100)','Training/Decision switch','Location','best');
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, sprintf('learning_curve_SNR%ddB.png',round(SNR_show))));
end

% constellations before/after
figure('Name','Constellations (one run)','Color','w');
subplot(1,2,1);
plot_constellation(y, S, 'Rx input y[n]');
subplot(1,2,2);
plot_constellation(xhat, S, 'Equalizer output \hat{x}[n]');

if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, sprintf('constellations_one_run_SNR%ddB.png',round(SNR_show))));
end

% SER vs SNR (Monte‑Carlo)
SNRdB = cfg.SNRdB_grid(:).';
SER   = zeros(size(SNRdB));
for iS = 0:numel(SNRdB)-1
    snrdb = SNRdB(iS+1);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
            simulate_one_run(cfg, h, D, snrdb);
        ser_mc(it) = ser_val;
        % ser_val
    end
    SER(iS+1) = mean(ser_mc);
    fprintf('SNR = %2d dB --> SER = %.3e\n', snrdb, SER(iS+1));
end

figure('Name','SER vs SNR','Color','w');
semilogy(SNRdB, SER, '-o','LineWidth',1.2); grid on;
xlabel('SNR (dB)'); ylabel('SER');
title(sprintf('SER vs SNR (N=%d, \\mu=%.4f, channel: %s)', ...
      cfg.equalizerLenN, cfg.mu, hDescr));
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_SNR.png'));
end


% sweep of step size mu
muGrid = [0.0005 0.001 0.002 0.003 0.005 0.01];
SER_mu = zeros(size(muGrid));
for im = 1:numel(muGrid)
    cfg2 = cfg; cfg2.mu = muGrid(im);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
            simulate_one_run(cfg2, h, D, cfg.plotOneRunSNRdB);
        ser_mc(it) = ser_val;
    end
    SER_mu(im) = mean(ser_mc);
    fprintf('mu = %.5f --> SER@%ddB = %.3e\n', cfg2.mu, cfg.plotOneRunSNRdB, SER_mu(im));
end

figure('Name','SER vs mu (at fixed SNR)','Color','w');
semilogy(muGrid, SER_mu, '-o','LineWidth',1.2); grid on;
xlabel('\mu'); ylabel('SER'); 
title(sprintf('SER vs \\mu at SNR=%.1f dB (N=%d)', cfg.plotOneRunSNRdB, cfg.equalizerLenN));
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_mu.png'));
end


% sweep of equalizer length N
Ngrid = [7 9 11 13 15];
SER_N = zeros(size(Ngrid));
for iN = 1:numel(Ngrid)
    cfg2 = cfg; cfg2.equalizerLenN = Ngrid(iN);
    D2 = pick_decision_delay(h, cfg2.equalizerLenN);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
            simulate_one_run(cfg2, h, D2, cfg.plotOneRunSNRdB);
        ser_mc(it) = ser_val;
    end
    SER_N(iN) = mean(ser_mc);
    fprintf('N = %d --> SER@%ddB = %.3e\n', cfg2.equalizerLenN, cfg.plotOneRunSNRdB, SER_N(iN));
end

figure('Name','SER vs equalizer length N','Color','w');
semilogy(Ngrid, SER_N, '-o','LineWidth',1.2); grid on;
xlabel('Equalizer length N'); ylabel('SER'); 
title(sprintf('SER vs N at SNR=%.1f dB (\\mu=%.4f)', cfg.plotOneRunSNRdB, cfg.mu));
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_N.png'));
end

%% benchmark test
% comprehensive joint parameter search (mu, N)
fprintf('\n=== Joint Grid Search: Optimal (mu, N) ===\n');
muGrid_full = [0.0005 0.001 0.002 0.003 0.005 0.007 0.01 0.015];
Ngrid_full = [5 7 9 11 13 15 17];
SNR_opt = cfg.plotOneRunSNRdB; % test at this SNR

SER_grid = zeros(numel(muGrid_full), numel(Ngrid_full));
fprintf('Testing %d mu values x %d N values = %d combinations...\n', ...
    numel(muGrid_full), numel(Ngrid_full), numel(muGrid_full)*numel(Ngrid_full));

for imu = 1:numel(muGrid_full)
    for iN = 1:numel(Ngrid_full)
        cfg_test = cfg;
        cfg_test.mu = muGrid_full(imu);
        cfg_test.equalizerLenN = Ngrid_full(iN);
        D_test = pick_decision_delay(h, cfg_test.equalizerLenN);
        
        ser_mc = zeros(1, cfg.numMC);
        for it = 1:cfg.numMC
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
                simulate_one_run(cfg_test, h, D_test, SNR_opt);
            ser_mc(it) = ser_val;
        end
        SER_grid(imu, iN) = mean(ser_mc);
        fprintf('  mu=%.5f, N=%2d --> SER = %.4e\n', ...
            cfg_test.mu, cfg_test.equalizerLenN, SER_grid(imu, iN));
    end
end

% find optimal combination
[minSER, idx] = min(SER_grid(:));
[opt_imu, opt_iN] = ind2sub(size(SER_grid), idx);
opt_mu = muGrid_full(opt_imu);
opt_N = Ngrid_full(opt_iN);

fprintf('\n*** OPTIMAL PARAMETERS ***\n');
fprintf('  mu = %.5f\n', opt_mu);
fprintf('  N  = %d\n', opt_N);
fprintf('  SER @ %ddB = %.4e\n', SNR_opt, minSER);
fprintf('*******************\n');

% visualize as heatmap
figure('Name','Joint Parameter Search','Color','w');
imagesc(Ngrid_full, 1:numel(muGrid_full), log10(SER_grid));
colorbar; colormap('hot');
set(gca, 'YTick', 1:numel(muGrid_full), 'YTickLabel', arrayfun(@(x) sprintf('%.5f',x), muGrid_full, 'UniformOutput', false));
xlabel('Equalizer Length N'); 
ylabel('\mu (step size)');
title(sprintf('log_{10}(SER) at SNR=%ddB - Joint Grid Search', SNR_opt));
hold on;
plot(opt_N, opt_imu, 'c*', 'MarkerSize', 20, 'LineWidth', 2);
text(opt_N, opt_imu, sprintf('  Optimal\n  (%.5f, %d)', opt_mu, opt_N), ...
    'Color', 'cyan', 'FontWeight', 'bold', 'FontSize', 10);
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'joint_grid_search.png'));
end

% plot SER vs N for each mu
figure('Name','SER vs N for different mu','Color','w');
for imu = 1:numel(muGrid_full)
    semilogy(Ngrid_full, SER_grid(imu,:), '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\mu=%.5f', muGrid_full(imu)));
    hold on;
end
grid on; xlabel('Equalizer Length N'); ylabel('SER');
title(sprintf('SER vs N for various \\mu at SNR=%ddB', SNR_opt));
legend('Location', 'best');
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_N_all_mu.png'));
end

% plot SER vs mu for each N
figure('Name','SER vs mu for different N','Color','w');
for iN = 1:numel(Ngrid_full)
    semilogy(muGrid_full, SER_grid(:,iN), '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('N=%d', Ngrid_full(iN)));
    hold on;
end
grid on; xlabel('\mu (step size)'); ylabel('SER');
title(sprintf('SER vs \\mu for various N at SNR=%ddB', SNR_opt));
legend('Location', 'best');
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_mu_all_N.png'));
end

% save workspace with results
save(fullfile(resultsDir, 'results.mat'));
fprintf('\nDone. Results saved under %s\n', fullfile('..','results'));
