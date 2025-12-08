% ECSE 512 - Digital Signal Processing (McGill)
% Term Project: Adaptive Equalization for 4-QAM over Dispersive Channel
clear; clc; close all;

cfg = struct();
% Simulation parameters
cfg.M = 4; % QAM order
cfg.sigma_s2 = 1;
cfg.trainLen = 1500;
cfg.dataLen = 15000;
cfg.equalizerLenN = 11; % equalizer length
cfg.mu = 0.005;
cfg.numMC = 5;
cfg.SNRdB_grid = 5:5:30;
cfg.saveFigs = true;
cfg.randomSeed = 42;
cfg.channel.type = 'severe';
cfg.channel.K = 5;
cfg.channel.randSTD = 1;
cfg.plotOneRunSNRdB = 30;
cfg.centerTap = [];
rng(cfg.randomSeed);

resultsDir = fullfile('..','results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

[S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));

% Channel setup
[h, hDescr] = generate_channel(cfg.channel);
h = normalize_channel(h);
K = numel(h);
fprintf('Channel (%s): h = %s\n', hDescr, mat2str(h,4));
if isempty(cfg.centerTap)
    D = pick_decision_delay(h, cfg.equalizerLenN);
else
    D = min(max(cfg.centerTap,0), cfg.equalizerLenN-1);
end
fprintf('Decision delay D = %d (equalizer length N = %d)\n', D, cfg.equalizerLenN);

% One run for learning curve
SNR_show = cfg.plotOneRunSNRdB;
fprintf('\n=== Single run at SNR = %.1f dB ===\n', SNR_show);
[x, y, w, noiseVar, xhat, e, d_used, c_hist, dec, mse, idxTrainEnd] = ...
    simulate_one_run(cfg, h, D, SNR_show);

figure('Name','Learning curve','Color','w');
plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth',1.2);
hold on;
yline(10*log10(mean(mse(end-100:end))), '--');
xline(idxTrainEnd,'--');
xlabel('n'); ylabel('E[|e[n]|^2] (dB)'); grid on;
title(sprintf('Learning curve, SNR=%.1f dB, N=%d, \\mu=%.4f, D=%d', ...
      SNR_show, cfg.equalizerLenN, cfg.mu, D));
legend('MSE (smoothed)','Steady-state','Training/Decision switch','Location','best');
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, sprintf('learning_curve_SNR%ddB.png',round(SNR_show))));
end

figure('Name','Constellations','Color','w');
subplot(1,2,1);
plot_constellation(y, S, 'Rx input y[n]');
subplot(1,2,2);
plot_constellation(xhat, S, 'Equalizer output \hat{x}[n]');

if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, sprintf('constellations_one_run_SNR%ddB.png',round(SNR_show))));
end

% SER vs SNR
SNRdB = cfg.SNRdB_grid(:).';
SER = zeros(size(SNRdB));
for iS = 1:numel(SNRdB)
    snrdb = SNRdB(iS);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = simulate_one_run(cfg, h, D, snrdb);
        ser_mc(it) = ser_val;
    end
    SER(iS) = mean(ser_mc);
    fprintf('SNR = %2d dB --> SER = %.3e\n', snrdb, SER(iS));
end

figure('Name','SER vs SNR','Color','w');
semilogy(SNRdB, SER, '-o','LineWidth',1.2); grid on;
xlabel('SNR (dB)'); ylabel('SER');
title(sprintf('SER vs SNR (N=%d, \\mu=%.4f, channel: %s)', ...
      cfg.equalizerLenN, cfg.mu, hDescr));
if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, 'SER_vs_SNR.png'));
end

% Sweep mu
muGrid = [0.0005 0.001 0.002 0.003 0.005 0.01];
SER_mu = zeros(size(muGrid));
for im = 1:numel(muGrid)
    cfg2 = cfg;
    cfg2.mu = muGrid(im);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = simulate_one_run(cfg2, h, D, cfg.plotOneRunSNRdB);
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

% Sweep N
Ngrid = [7 9 11 13 15];
SER_N = zeros(size(Ngrid));
for iN = 1:numel(Ngrid)
    cfg2 = cfg;
    cfg2.equalizerLenN = Ngrid(iN);
    D2 = pick_decision_delay(h, cfg2.equalizerLenN);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = simulate_one_run(cfg2, h, D2, cfg.plotOneRunSNRdB);
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

save(fullfile(resultsDir, 'results.mat'));
fprintf('\nDone. Results saved under %s\n', fullfile('..','results'));
