function [xhat, e, c_hist, d_used, decisions, mse, idxTrainEnd] = lms_equalizer(y, x, N, mu, Ltrain, D, S)
% LMS adaptive equalizer with training then decision-directed mode
    y = y(:);
    x = x(:);
    L = numel(y);
    S = S(:).';
    xhat = complex(zeros(L,1));
    e = complex(zeros(L,1));
    decisions = complex(zeros(L,1));
    d_used = complex(zeros(L,1));
    mse = zeros(L,1);
    c = zeros(N,1);
    c(min(D+1,N)) = 1; % initialize at delay D
    c_hist = complex(zeros(N,L));

    idxStart = N;
    idxTrainEnd = min(L, idxStart + Ltrain - 1);

    for n = idxStart:L
        yvec = flipud(y(n-N+1:n));
        xhat(n) = c' * yvec;
        if n <= idxTrainEnd
            d = 0;
            if (n-D) >= 1 && (n-D) <= numel(x)
                d = x(n-D);
            end
        else
            [~, d] = qam4_demod(xhat(n), sqrt(mean(abs(S).^2)));
        end
        d_used(n) = d;
        e(n) = d - xhat(n);
        c = c + mu * yvec * conj(e(n));
        c_hist(:,n) = c;
        decisions(n) = d;
        mse(n) = abs(e(n)).^2;
    end
end
