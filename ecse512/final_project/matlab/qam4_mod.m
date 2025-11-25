function x = qam4_mod(idx, sigma_s)
%QAM4_MOD Maps integer indices 1..4 to 4‑QAM symbols.
%   x = qam4_mod(idx, sigma_s)
%     idx       vector with entries in {1,2,3,4}
%     sigma_s   amplitude scaling (so E{|x|^2} = sigma_s^2)
    if nargin < 2 || isempty(sigma_s), sigma_s = 1; end
    [S,~] = qam4_constellation(sqrt(sigma_s^2));
    x = S(idx);
end
