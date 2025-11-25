function [y_noisy, w, sigma2] = add_awgn(y_clean, SNRdB, signalPower)
%ADD_AWGN Adds complex circular AWGN to reach the target SNR at the input.
%   signalPower is E{|signal|^2}. For our system, E{|y_clean|^2} ~= signalPower
%   because the channel is unit energy and input symbol power is sigma_s2.
    if nargin < 3 || isempty(signalPower)
        signalPower = mean(abs(y_clean).^2);
    end
    SNRlin = 10.^(SNRdB/10);
    sigma2 = signalPower ./ SNRlin;
    w = sqrt(sigma2/2) * (randn(size(y_clean)) + 1j*randn(size(y_clean)));
    y_noisy = y_clean + w;
end
