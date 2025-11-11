
%% ECSE 563 — Assignment 3 — Build & Run All Parts
% Open MATLAB in this folder, run:
%   A3_params; A3_run_all

A3_params;  % load parameters

%% ---------- Part 1: primary control only ----------
mdl1 = A3_part1_build();
set_param([mdl1 '/RateLimit1'], 'Commented', 'on'); % ensure disabled
set_param(mdl1, 'StopTime', '30');
out1 = sim(mdl1);

% Save frequency and mechanical powers
f1 = out1.f_Hz; Pm1 = out1.Pm_MW;

%% ---------- Part 1(e): ramp limit on Gen 1 ----------
set_param([mdl1 '/RateLimit1'], 'Commented', 'off'); % enable 1 MW/s
set_param(mdl1, 'StopTime', '120');                 % give it more time
out1e = sim(mdl1);
fe = out1e.f_Hz; Pm1e = out1e.Pm_MW;

%% ---------- Part 2: add secondary control (AGC) ----------
mdl2 = A3_part2_build();
set_param([mdl2 '/RateLimit1'], 'Commented', 'on');  % no ramp limit
set_param(mdl2, 'StopTime', '900');
out2 = sim(mdl2);
f2 = out2.f_Hz; Pm2 = out2.Pm_MW;

%% ---------- Part 3: low inertia / damping ----------
M = 0.6*M; D = 0.5*D;   % adjust parameters
mdl3 = A3_part3_build();
set_param([mdl3 '/RateLimit1'], 'Commented', 'on');
set_param(mdl3, 'StopTime', '900');
out3 = sim(mdl3);
f3 = out3.f_Hz; Pm3 = out3.Pm_MW;

%% ---------- Save plots ----------
% Build helper to get timeseries data
tt = f1.time; f1sig = f1.signals.values;
tPm = Pm1.time; Pm1sig = Pm1.signals.values;  % Nx2 matrix

% Part 1 plots
figure; plot(f1.time, f1.signals.values); grid on; xlabel('Time [s]'); ylabel('\Delta f [Hz]');
title('Part 1: Frequency deviation'); saveas(gcf, 'part1_freq.png');

figure; plot(Pm1.time, Pm1.signals.values(:,1),'LineWidth',1); hold on;
plot(Pm1.time, Pm1.signals.values(:,2),'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('\Delta P_m [MW]'); legend('Gen 1','Gen 2','Location','best');
title('Part 1: Mechanical powers'); saveas(gcf, 'part1_powers.png');

% Part 1(e)
figure; plot(fe.time, fe.signals.values, 'LineWidth',1); grid on; xlabel('Time [s]'); ylabel('\Delta f [Hz]');
title('Part 1(e): Frequency with 1 MW/s limit on Gen 1'); saveas(gcf, 'part1e_freq.png');

figure; plot(Pm1e.time, Pm1e.signals.values(:,1),'LineWidth',1); hold on;
plot(Pm1e.time, Pm1e.signals.values(:,2),'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('\Delta P_m [MW]'); legend('Gen 1','Gen 2'); 
title('Part 1(e): Mechanical powers'); saveas(gcf, 'part1e_powers.png');

% Part 2
figure; plot(f2.time, f2.signals.values,'LineWidth',1); grid on; xlabel('Time [s]'); ylabel('\Delta f [Hz]');
title('Part 2: Frequency with AGC'); saveas(gcf, 'part2_freq.png');

% Part 3
figure; plot(f3.time, f3.signals.values,'LineWidth',1); grid on; xlabel('Time [s]'); ylabel('\Delta f [Hz]');
title('Part 3: Frequency with lower M and D'); saveas(gcf, 'part3_freq.png');

disp('All simulations complete. Figures saved in the current folder.');
