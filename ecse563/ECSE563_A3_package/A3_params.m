
%% ECSE 563 — Assignment 3 — Parameters
%% Units: MW, Hz, seconds
f0     = 60;             % nominal frequency [Hz]
Psys   = 2000;           % system base [MW]
DeltaP = 100;         % load step [MW]
D      = 20;              % load frequency sensitivity [MW/Hz]
M      = 700;              % inertia parameter in MW*s/Hz (M = 2*E_k/f0)

% Governor & turbine time constants
TG1 = 0.25;  TCH1 = 3;
TG2 = 0.5;  TCH2 = 8;

% Droop (speed regulation) in Hz/MW
R1 = 0.006;     % 3 Hz / 500 MW
R2 = 0.012;     % 3 Hz / 250 MW
invR1 = 1/R1; invR2 = 1/R2;   % MW/Hz

% Secondary control (AGC)
Ki_AGC = 1;      % MW/(Hz*s), tuned to ~10–15 min recovery
a1 = 0.75; a2 = 0.25; % participation factors (sum to 1)
