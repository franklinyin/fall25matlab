function [idx, decisions] = qam4_demod(z, sigma_s)
% Minimum‑distance decision device for 4‑QAM.
    if nargin < 2 || isempty(sigma_s), sigma_s = 1; end
    [S,~] = qam4_constellation(sqrt(sigma_s^2));
    idx = zeros(size(z));
    decisions = zeros(size(z));
    for n = 1:numel(z)
        [~,k] = min(abs(z(n) - S));
        idx(n) = k;
        decisions(n) = S(k);
    end
end
