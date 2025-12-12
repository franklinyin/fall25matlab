% main script

clear; clc; close all;

cfg = struct();

% simulation parameters
cfg.M = 4;
cfg.sigma_s2 = 1;
cfg.trainLen = 1500;
cfg.dataLen = 15000;
cfg.equalizerLenN = 9;
cfg.mu = 0.015;
cfg.randomSeed = 42;
cfg.saveFigs = true;

rng(cfg.randomSeed);

% make output dir
resultsDir = fullfile('..','results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

%% constellation diagrams
fprintf('\nconstellations...\n');

channelTypes = {'low', 'mild', 'severe'};
SNRdBs = [10, 20, 30];

for chIdx = 1:length(channelTypes)
    channelType = channelTypes{chIdx};
    cfg.channel.type = channelType;
    [h, hDescr] = generate_channel(cfg.channel);
    h = normalize_channel(h);
    D = pick_decision_delay(h, cfg.equalizerLenN);
    fprintf('  %s\n', hDescr);
    for snrIdx = 1:length(SNRdBs)
        snr = SNRdBs(snrIdx);
        plot_constellation_comparison(cfg, h, hDescr, D, snr, resultsDir);
    end
end

%% learning curves: mild at different SNRs
fprintf('\nlearning curves (mild, SNR sweep)...\n');
plot_learning_curve_comparison(cfg, {'mild'}, [10, 20, 30], resultsDir, 'snr');

%% learning curves: channels at 20dB
fprintf('learning curves (channel sweep, 20dB)...\n');
plot_learning_curve_comparison(cfg, {'low', 'mild', 'severe'}, [20], resultsDir, 'channel');

%% learning curves: mu sweep
fprintf('learning curves (mu sweep)...\n');
muGrid = [0.001, 0.01, 0.05, 0.1];
plot_learning_curve_mu_comparison(cfg, 'mild', 20, muGrid, resultsDir);

%% learning curves: N sweep
fprintf('learning curves (N sweep)...\n');
Ngrid = [3, 7, 11, 17];
plot_learning_curve_N_comparison(cfg, 'mild', 20, Ngrid, resultsDir);

%% SER vs SNR
fprintf('SER vs SNR...\n');
cfg.numMC = 5;
SNRdB_grid = 0:3:30;
plot_ser_snr_comparison(cfg, {'low', 'mild', 'severe'}, SNRdB_grid, resultsDir);

fprintf('\ndone. saved to %s\n', resultsDir);