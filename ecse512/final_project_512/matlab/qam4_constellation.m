function [S, sigma_s] = qam4_constellation(sigma_s)
% Returns 4-QAM constellation
    if nargin < 1 || isempty(sigma_s)
        sigma_s = 1;
    end
    S = sigma_s/sqrt(2) * [ 1+1j, -1+1j, -1-1j, 1-1j ];
end
