function [x, y, w, noiseVar, xhat, e, d_used, c_hist, decisions, mse, idxTrainEnd, ser] = ...
    simulate_one_run(cfg, h, D, SNRdB)
% Generate one Monte‑Carlo realization and equalize.
%   Returns many internal signals for analysis/plotting.

    L = cfg.trainLen + cfg.dataLen;
    % IID 4‑QAM indices 1..4
    idx = randi(4, L, 1);
    x = qam4_mod(idx, sqrt(cfg.sigma_s2));

    % channel + noise (unit energy channel -> output power ~ sigma_s2)
    y_clean = filter(h, 1, x);
    % mean(abs(y_clean).^2)
    [y, w, noiseVar] = add_awgn(y_clean, SNRdB, cfg.sigma_s2);

    % equalization (training then decision‑directed)
    N = cfg.equalizerLenN;
    mu = cfg.mu;
    [xhat, e, c_hist, d_used, decisions, mse, idxTrainEnd] = ...
        lms_equalizer(y, x, N, mu, cfg.trainLen, D, qam4_constellation(sqrt(cfg.sigma_s2)));

    % SER over the decision‑directed phase only, compare to aligned truth x[n-D]
    n0 = idxTrainEnd + 1;
    if n0 < 1 || n0 > L
        ser = NaN;
    else
        x_true = zeros(L-n0+1,1);
        x_hard = zeros(L-n0+1,1);
        for n = n0:L
            if (n-D)>=1 && (n-D)<=L
                x_true(n-n0+1) = x(n-D);
            else
                x_true(n-n0+1) = 0;
            end
            x_hard(n-n0+1) = decisions(n);
        end
        % map to indices to avoid comparing complex floats
        [idx_true,~] = qam4_demod(x_true, sqrt(cfg.sigma_s2));
        [idx_hard,~] = qam4_demod(x_hard, sqrt(cfg.sigma_s2));
        % sum(idx_true~=idx_hard)
        ser = measure_ser(idx_true, idx_hard);
    end
end
