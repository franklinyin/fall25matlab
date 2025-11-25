% parameters
M = 48;
N = M + 1;
beta = 3.68;
n = 0:M;


hd = sin(0.3*pi*(n-24))./(pi*(n-24)) + ...
     0.5*(-1).^(n-24).*sin(0.4*pi*(n-24))./(pi*(n-24));
hd(n==24) = 0.3 + 0.5*0.4;

% window and FIR coefficients 
w  = kaiser(N, beta).';
h  = hd .* w;

% plot magnetude response---
[H, wgrid] = freqz(h, 1, 4096);
figure; plot(wgrid/pi, abs(H)); grid on
xlabel('\omega/\pi'); ylabel('|H(e^{j\omega})|');
title('Q7.44 part d (M=48, \beta=3.68)');

% Save plot to PNG
saveas(gcf, 'hw5_7_44_d.png');
