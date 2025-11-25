clc; clear; close all;

M = 48;
N = M + 1;
beta = 3.68;
B = 1.0;
C = 0.5;

w1 = 0.3*pi;
w2 = 0.6*pi;

A_dB = 42.4;
Dw = (A_dB - 8) / (2.285 * N);
wp1 = 0.3*pi - Dw/2; ws1 = 0.3*pi + Dw/2;
ws2 = 0.6*pi - Dw/2; wp2 = 0.6*pi + Dw/2;

%% result from part c
n = 0:M;
k = n - M/2; % centered index (M/2=24)
hd = sin(0.3*pi*(n-24))./(pi*(n-24)) + ...
     0.5*(-1).^(n-24).*sin(0.4*pi*(n-24))./(pi*(n-24));
hd(n==24) = 0.3 + 0.5*0.4;

%% part d, choosing Kaiser window design
wK = kaiser(N, beta).';
hK = hd .* wK;

[Hk, wgrid] = freqz(hK, 1, 32768);
magK = abs(Hk);
magK_dB = 20*log10(max(magK, 1e-12));


% Linear magnitude
figure('Name','Linear Magnitude');
plot(wgrid/pi, magK, 'LineWidth',1.2); hold on;
grid on; xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Kaiser Magnitude Response (linear magnitude)');
legend('Kaiser-window', 'Parks–McClellan (firpm)', 'Location','Best');
% Save plot to PNG
saveas(gcf, 'hw5_7_44_d_linear.png');

% dB magnitude
figure('Name','Magnitude in dB');
plot(wgrid/pi, magK_dB, 'LineWidth',1.2); hold on;
grid on; xlabel('\omega/\pi'); ylabel('Magnitude (dB)');
title('Kaiser Magnitude Response (dB)');
legend('Kaiser-window', 'Parks–McClellan (firpm)', 'Location','Best');
% Save plot to PNG
saveas(gcf, 'hw5_7_44_d_db.png');

%% part e getting numerical Measure δ1, δ2, δ3 on spec bands
% Spec bands from wp1, ws1, ws2, wp2
idx_p1 = find(wgrid >= 0 & wgrid <= wp1);
idx_s = find(wgrid >= ws1 & wgrid <= ws2);
idx_p2 = find(wgrid >= wp2 & wgrid <= pi);

delta1_K = max(abs(magK(idx_p1) - B));
delta2_K = max(magK(idx_s));
delta3_K = max(abs(magK(idx_p2) - C));

fprintf('part e');
fprintf('  delta1 (low passband) = %.6g\n', delta1_K);
fprintf('  delta2 (stopband) = %.6g\n', delta2_K);
fprintf('  delta3 (high passband) = %.6g\n', delta3_K);
fprintf('  (Approx attenuation = %.1f dB → delta ≈ %.3g)\n', A_dB, 10^(-A_dB/20));
fprintf('\n');

%% part f, choosing firpm (Parks–McClellan)
% Frequency grid for firpm is normalized to 1 ↔ π (use /pi).
f = [0  wp1 ws1 ws2 wp2 pi]/pi;
a = [B B 0 0 C C];

% Weight selection: Balance relative ripple (weight high band ~ 2x)
wts = [1 1 2];  % [pass1, stop, pass2]

% Equiripple design
hPM = firpm(M, f, a, wts);

[Hpm, wgrid2] = freqz(hPM, 1, 32768);
magPM = abs(Hpm);
magPM_dB = 20*log10(max(magPM, 1e-12));

% Measure firpm ripples on the same spec bands
delta1_PM = max(abs(magPM(idx_p1) - B));
delta2_PM = max(magPM(idx_s));
delta3_PM = max(abs(magPM(idx_p2) - C));

%% ---------------- (g) Plots & comparison ----------------
% Linear magnitude
figure('Name','Linear Magnitude');
plot(wgrid/pi, magK, 'LineWidth',1.2); hold on;
plot(wgrid2/pi, magPM, 'LineWidth',1.2);
grid on; xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Kaiser vs. Parks–McClellan (linear magnitude)  [firpm weights = [1 1 2]]');
legend('Kaiser-window', 'Parks–McClellan (firpm)', 'Location','Best');
% Save plot to PNG
saveas(gcf, 'hw5_7_44_f_linear.png');

% dB magnitude
figure('Name','Magnitude in dB');
plot(wgrid/pi, magK_dB, 'LineWidth',1.2); hold on;
plot(wgrid2/pi, magPM_dB, 'LineWidth',1.2);
grid on; xlabel('\omega/\pi'); ylabel('Magnitude (dB)');
title('Kaiser vs. Parks–McClellan (dB)  [firpm weights = [1 1 2]]');
legend('Kaiser-window', 'Parks–McClellan (firpm)', 'Location','Best');

% Annotate spec bands on the dB plot (optional helper lines)
yl = ylim;
for xx = [wp1, ws1, ws2, wp2]/pi
    xline(xx, ':', 'Color',[0.5 0.5 0.5]);
end
ylim(yl);

% Save plot to PNG
saveas(gcf, 'hw5_7_44_f_db.png');

%% ---------------- Console summary ----------------
fprintf('=== SPEC BANDS (from Kaiser Δω estimate) ===\n');
fprintf('  wp1/pi = %.6f\n', wp1/pi);
fprintf('  ws1/pi = %.6f\n', ws1/pi);
fprintf('  ws2/pi = %.6f\n', ws2/pi);
fprintf('  wp2/pi = %.6f\n', wp2/pi);
fprintf('\n');

fprintf('=== (e) Kaiser-window measured ripples ===\n');
fprintf('  delta1 (low passband) = %.6g\n', delta1_K);
fprintf('  delta2 (stopband) = %.6g\n', delta2_K);
fprintf('  delta3 (high passband) = %.6g\n', delta3_K);
fprintf('  (Approx attenuation = %.1f dB → delta ≈ %.3g)\n', A_dB, 10^(-A_dB/20));
fprintf('\n');

fprintf('=== (f) firpm measured ripples (weights = [%g %g %g]) ===\n', wts(1), wts(2), wts(3));
fprintf('  delta1 (low passband) = %.6g\n', delta1_PM);
fprintf('  delta2 (stopband) = %.6g\n', delta2_PM);
fprintf('  delta3 (high passband) = %.6g\n', delta3_PM);
fprintf('\n');

% High-level comparison
fprintf('=== (g) Brief comparison ===\n');
fprintf('- Kaiser-window: quick, predictable Δω & δ from (beta,N); not equiripple.\n');
fprintf('- firpm: equiripple; generally achieves smaller max error per band for same M.\n');
fprintf('- Using weights [1 1 2] balances relative passband ripples (since high band level is 0.5).\n');
