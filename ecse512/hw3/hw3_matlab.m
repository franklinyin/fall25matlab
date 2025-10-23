% matlab plot problem
a = 0.8;              % choose |a|<1 for convergence
N_values = [16, 32, 64];

for N = N_values
    k = 0:N-1;
    omega_k = 2*pi*k/N;
    Xk = 1 ./ (1 - a*exp(-1j*omega_k));  % sampled DTFT derived in the pdf
    x_hat = ifft(Xk);                    % IDFT
    
    figure;
    stem(0:N-1, real(x_hat), 'filled');
    title(['Reconstructed x[n] from N = ', num2str(N)]);
    xlabel('n');
    ylabel('Re{x[n]}');
    
    % saving as imgs
    filename = sprintf('Reconstructed_x_N%d.png', N);
    saveas(gcf, filename);
    
end
