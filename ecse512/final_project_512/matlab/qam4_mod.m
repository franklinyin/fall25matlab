function x = qam4_mod(idx, sigma_s)
% Maps indices 1..4 to 4-QAM symbols
    if nargin < 2 || isempty(sigma_s)
        sigma_s = 1;
    end
    [S,~] = qam4_constellation(sqrt(sigma_s^2));
    x = S(idx);
end
