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

%% Generate constellation diagrams for all combinations
fprintf('=== Generating Constellation Diagrams ===\n');

channelTypes = {'low', 'mild', 'severe'};
SNRdBs = [10, 20, 30];

for chIdx = 1:length(channelTypes)
    channelType = channelTypes{chIdx};
    cfg.channel.type = channelType;
    
    % Generate and normalize channel
    [h, hDescr] = generate_channel(cfg.channel);
    h = normalize_channel(h);
    D = pick_decision_delay(h, cfg.equalizerLenN);
    
    fprintf('\nChannel: %s\n', hDescr);
    
    for snrIdx = 1:length(SNRdBs)
        snr = SNRdBs(snrIdx);
        fprintf('  Generating constellation for SNR = %d dB... ', snr);
        
        plot_constellation_comparison(cfg, h, hDescr, D, snr, resultsDir);
    end
end

%% Generate learning curve comparison: mild channel at different SNRs
fprintf('\n=== Generating Learning Curve: Mild Channel at 10/20/30 dB ===\n');

plot_learning_curve_comparison(cfg, {'mild'}, [10, 20, 30], resultsDir, 'snr');

%% Generate learning curve comparison: different channels at 20 dB
fprintf('\n=== Generating Learning Curve: Low/Mild/Severe at 20 dB ===\n');

plot_learning_curve_comparison(cfg, {'low', 'mild', 'severe'}, [20], resultsDir, 'channel');


%% Generate learning curve comparison for different step sizes
fprintf('\n=== Generating Learning Curve: Different Step Sizes ===\n');

muGrid = [0.001, 0.01, 0.05, 0.1];
plot_learning_curve_mu_comparison(cfg, 'mild', 20, muGrid, resultsDir);

%% Generate learning curve comparison for different equalizer lengths
fprintf('\n=== Generating Learning Curve: Different Equalizer Lengths ===\n');

Ngrid = [3, 7, 11, 17];
plot_learning_curve_N_comparison(cfg, 'mild', 20, Ngrid, resultsDir);


%% Generate SER vs SNR comparison for all channels
fprintf('\n=== Generating SER vs SNR Comparison: Low/Mild/Severe ===\n');

cfg.numMC = 5; % Monte Carlo trials for SER
SNRdB_grid = 0:3:30; % SNR range for comparison

plot_ser_snr_comparison(cfg, {'low', 'mild', 'severe'}, SNRdB_grid, resultsDir);

fprintf('\n=== All plots generated successfully! ===\n');
fprintf('Results saved in: %s\n', resultsDir);