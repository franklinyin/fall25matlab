function [y_noisy, w, sigma2] = add_awgn(y_clean, SNRdB, signalPower)
% Adds complex AWGN to reach target SNR
    if nargin < 3 || isempty(signalPower)
        signalPower = mean(abs(y_clean).^2);
    end
    SNRlin = 10.^(SNRdB/10);
    sigma2 = signalPower ./ SNRlin;
    w = sqrt(sigma2/2) * (randn(size(y_clean)) + 1j*randn(size(y_clean)));
    y_noisy = y_clean + w;
end
