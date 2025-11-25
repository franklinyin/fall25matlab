# ECSE 512 – Adaptive 4‑QAM Equalization (MATLAB)

This repository contains a complete MATLAB implementation of an adaptive FIR equalizer trained with the complex LMS algorithm for 4‑QAM transmission over an ISI channel, plus a LaTeX report template.

## How to run

1. Open MATLAB in the `matlab/` folder.
2. Run `main.m`. It will:
   - Generate 4‑QAM symbols (training + payload);
   - Filter them through a normalized ISI channel and add complex AWGN for a set of SNRs;
   - Train an LMS equalizer (training then decision‑directed);
   - Produce and save figures under `../results/` and a `.mat` file with all results.

You can change simulation parameters at the top of `main.m` (filter lengths, step size, SNR grid, training length, Monte‑Carlo runs, channel type, etc.).

## MATLAB version

The code uses only base MATLAB functions (`filter`, `conv`, `randn`, plotting), so no toolboxes are required.
