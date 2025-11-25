% ECSE 512 – Digital Signal Processing (McGill)
% Term Project: Adaptive Equalization for 4‑QAM over Dispersive Channel
% Main driver script
% -------------------------------------------------------------------------
clear; clc; close all;

cfg = struct();
% --------------------- Simulation parameters ------------------------------
cfg.M                = 4;           % QAM order (fixed to 4‑QAM)
cfg.sigma_s2         = 1;           % symbol power (E{|x|^2})
cfg.trainLen         = 1500;        % training symbols
cfg.dataLen          = 15000;       % data symbols (decision‑directed phase)
cfg.equalizerLenN    = 11;          % LMS equalizer length (N taps)
cfg.mu               = 0.003;       % LMS step size (tune per scenario)
cfg.numMC            = 5;           % Monte‑Carlo trials (increase for stats)
cfg.SNRdB_grid       = 5:5:30;      % SNR sweep for SER curves
cfg.saveFigs         = true;        % save figures to ../results
cfg.randomSeed       = 42;          % reproducibility
cfg.channel.type     = 'severe';    % 'mild' | 'severe' | 'random'
cfg.channel.K        = 5;           % channel length for 'random' type
cfg.channel.randSTD  = 1;           % std for Rayleigh taps before normalization
cfg.plotOneRunSNRdB  = 20;          % SNR for which to show a detailed run
cfg.centerTap        = [];          % optional target decision delay D (0..N-1)
% -------------------------------------------------------------------------
rng(cfg.randomSeed);

% Make output dir
resultsDir = fullfile('..','results');
if ~exist(resultsDir,'dir'); mkdir(resultsDir); end

% Reference 4‑QAM constellation (unit power scaled to cfg.sigma_s2)
[S, ~] = qam4_constellation(sqrt(cfg.sigma_s2));  % S is 4 symbols

% ------------------- Channel definition and decision delay ----------------
[h, hDescr] = generate_channel(cfg.channel);
h = normalize_channel(h);           % unit energy
K = numel(h);
fprintf('Channel (%s): h = %s\n', hDescr, mat2str(h,4));

% Choose decision delay D to align the equalizer target
if isempty(cfg.centerTap)
    D = pick_decision_delay(h, cfg.equalizerLenN);
else
    D = min(max(cfg.centerTap,0), cfg.equalizerLenN-1);
end
fprintf('Decision delay D = %d (equalizer length N = %d)\n', D, cfg.equalizerLenN);

% ------------------- One illustrative run (learning curve) ----------------
SNR_show = cfg.plotOneRunSNRdB;
fprintf('\n=== Single run at SNR = %.1f dB for learning curve ===\n', SNR_show);
[x, y, w, noiseVar, xhat, e, d_used, c_hist, dec, mse, idxTrainEnd] = ...
    simulate_one_run(cfg, h, D, SNR_show);

% Learning curve (averaged MSE over Monte‑Carlo is plotted later; here just one run)
figure('Name','Learning curve (one run)','Color','w');
plot(1:numel(mse), 10*log10(movmean(mse, 50)), 'LineWidth',1.2);
hold on; yline(10*log10(mean(mse(end-100:end))), '--');
xline(idxTrainEnd,'--');
xlabel('n'); ylabel('E[|e[n]|^2] (dB)'); grid on;
title(sprintf('Learning curve, SNR=%.1f dB, N=%d, \\mu=%.4f, D=%d', ...
      SNR_show, cfg.equalizerLenN, cfg.mu, D));
legend('MSE (smoothed)','Steady‑state (avg last 100)','Training/Decision switch','Location','best');

% Constellations before/after
figure('Name','Constellations (one run)','Color','w');
subplot(1,2,1);
plot_constellation(y, S, 'Rx input y[n]');
subplot(1,2,2);
plot_constellation(xhat, S, 'Equalizer output \hat{x}[n]');

if cfg.saveFigs
    saveas(gcf, fullfile(resultsDir, sprintf('constellations_one_run_SNR%ddB.png',round(SNR_show))));
end

% ------------------- SER vs SNR (Monte‑Carlo) -----------------------------
SNRdB = cfg.SNRdB_grid(:).';
SER   = zeros(size(SNRdB));
for iS = 0:numel(SNRdB)-1
    snrdb = SNRdB(iS+1);
    ser_mc = zeros(1, cfg.numMC);
    for it = 1:cfg.numMC
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ser_val] = ...
            simulate_one_run(cfg, h, D, snrdb);
        ser_mc(it) = ser_val;
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

% ------------------- Sweep of step size mu (optional) ---------------------
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

% ------------------- Sweep of equalizer length N (optional) ---------------
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

% Save workspace with results
save(fullfile(resultsDir, 'results.mat'));
fprintf('\nDone. Results saved under %s\n', fullfile('..','results'));
