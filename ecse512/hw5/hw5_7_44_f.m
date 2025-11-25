%  Bands (normalized to pi) 
A   = 42.423;                        % from beta=3.68
Dw  = (A - 8) / (2.285*(M+1));      % transition width (rad)
wp1 = 0.3*pi - Dw/2;  ws1 = 0.3*pi + Dw/2;
ws2 = 0.6*pi - Dw/2;  wp2 = 0.6*pi + Dw/2;

f   = [0  wp1 ws1 ws2 wp2 pi]/pi;   % edges (0..1)
a   = [1  1   0   0   0.5 0.5];     % desired levels

% Use weights inversely proportional to the ripples from (e)
d1 = 0.01095; d2 = 0.01081; d3 = 0.00495;
wts = [1/d1  1/d2  1/d3];

%  Equiripple design 
h_pm = firpm(M, f, a, wts);

%  Plot 
[Hpm, wgrid] = freqz(h_pm, 1, 4096);
figure; plot(wgrid/pi, abs(Hpm)); grid on
xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Parks–McClellan multiband FIR (firpm, M=48)');

% Save plot to PNG
saveas(gcf, 'hw5_7_44_f.png');
