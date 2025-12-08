function [x, y, w, noiseVar, xhat, e, d_used, c_hist, decisions, mse, idxTrainEnd, ser] = simulate_one_run(cfg, h, D, SNRdB)
% One simulation run
    L = cfg.trainLen + cfg.dataLen;
    idx = randi(4, L, 1);
    x = qam4_mod(idx, sqrt(cfg.sigma_s2));

    y_clean = filter(h, 1, x);
    [y, w, noiseVar] = add_awgn(y_clean, SNRdB, cfg.sigma_s2);

    N = cfg.equalizerLenN;
    mu = cfg.mu;
    [xhat, e, c_hist, d_used, decisions, mse, idxTrainEnd] = lms_equalizer(y, x, N, mu, cfg.trainLen, D, qam4_constellation(sqrt(cfg.sigma_s2)));

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
        [idx_true,~] = qam4_demod(x_true, sqrt(cfg.sigma_s2));
        [idx_hard,~] = qam4_demod(x_hard, sqrt(cfg.sigma_s2));
        ser = measure_ser(idx_true, idx_hard);
    end
end
