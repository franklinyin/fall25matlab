function [xhat, e, c_hist, d_used, decisions, mse, idxTrainEnd] = lms_equalizer(y, x, N, mu, Ltrain, D, S)
%LMS_EQUALIZER Complex LMS adaptive FIR equalizer with training then DD.
%   Implements:
%       xhat[n] = sum_{k=0}^{N-1} c_k[n] * y[n-k]
%       e[n]    = d[n] - xhat[n]
%       c[n+1]  = c[n] + mu * y_vec[n] * conj(e[n])
%   where y_vec[n] = [y[n], y[n-1], ..., y[n-N+1]]^T
%
%   Inputs:
%       y      : received samples (1xL or Lx1)
%       x      : transmitted symbol sequence aligned as d[n] = x[n-D]
%       N      : number of equalizer taps
%       mu     : LMS step size
%       Ltrain : number of training iterations
%       D      : decision delay (0..N-1)
%       S      : 4‑QAM constellation (1x4)
%
%   Outputs:
%       xhat        : equalized output (1xL)
%       e           : error signal used by LMS (1xL)
%       c_hist      : N x L evolution of the equalizer taps
%       d_used      : desired sequence actually used at each n
%       decisions   : hard decisions in DD mode
%       mse         : |e[n]|^2 (vector)
%       idxTrainEnd : index where training stopped (N-1 + Ltrain clipped)
%
%   Notes:
%       - Updates start at n = N (1‑based) so that y[n-k] is defined.
%       - During training, d[n] = x[n-D].
%       - During decision‑directed mode, d[n] = Q( xhat[n] ).
%       - The first N-1 outputs are left as zero.

    y = y(:); x = x(:); L = numel(y);
    S = S(:).'; % row
    xhat = complex(zeros(L,1));
    e    = complex(zeros(L,1));
    decisions = complex(zeros(L,1));
    d_used = complex(zeros(L,1));
    mse  = zeros(L,1);
    c    = zeros(N,1);
    % initialize as delta at the chosen delay for faster convergence
    c(min(D+1,N)) = 1;
    c_hist = complex(zeros(N,L));

    idxStart = N; % first index where all y[n-k] are known
    idxTrainEnd = min(L, idxStart + Ltrain - 1);

    for n = idxStart:L
        yvec = flipud(y(n-N+1:n));        % [y[n], y[n-1], ..., y[n-N+1]]
        xhat(n) = c' * yvec;
        if n <= idxTrainEnd
            d = 0;
            if (n-D) >= 1 && (n-D) <= numel(x)
                d = x(n-D);
            end
        else
            % Decision‑directed
            [~, d] = qam4_demod(xhat(n), sqrt(mean(abs(S).^2)));
        end
        d_used(n) = d;
        e(n) = d - xhat(n);
        c    = c + mu * yvec * conj(e(n));
        c_hist(:,n) = c;
        decisions(n) = d;  % in training equals reference x (delayed)
        mse(n) = abs(e(n)).^2;
    end
end
