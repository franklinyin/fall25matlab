% Parameters
M = 64;  % window length - 1 (so total points = M+1)
n = 0:M;

% define Hanning window
w = 0.5 * (1 - cos(2*pi*n/M));

% compute FT
NFFT = 1024;  % High resolution FFT
W = fftshift(fft(w, NFFT));
omega = linspace(-pi, pi, NFFT);

% normalize magnitude (in dB)
W_mag = 20*log10(abs(W) / max(abs(W)));

% plot
figure;
plot(omega, W_mag, 'LineWidth', 1.5);
xlabel('\omega (radians/sample)');
ylabel('Magnitude (dB)');
title('Normalized Magnitude Spectrum of Hanning Window');
grid on;
axis tight;
xlim([0,pi]);
ylim([-100 0]);

% Save plot to PNG
saveas(gcf, 'hw5_7_43_c_sketch.png');
