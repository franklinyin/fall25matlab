function [S, sigma_s] = qam4_constellation(sigma_s)
%QAM4_CONSTELLATION Returns 4‑QAM (QPSK) constellation scaled by sigma_s.
%   S = sigma_s/sqrt(2) * [ 1+1j, -1+1j, -1-1j, 1-1j ]
%   S: 1x4 complex array with symbols
%   sigma_s: scalar amplitude scaling (defaults to 1)
    if nargin < 1 || isempty(sigma_s)
        sigma_s = 1;
    end
    S = sigma_s/sqrt(2) * [ 1+1j, -1+1j, -1-1j, 1-1j ];
end
