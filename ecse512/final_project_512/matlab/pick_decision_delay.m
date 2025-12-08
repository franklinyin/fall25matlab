function D = pick_decision_delay(h, N)
% Choose a reasonable decision delay 0..N-1.
%   Heuristic: place center tap near the largest‑energy channel tap.
    [~, idx] = max(abs(h));
    % D = idx - 1;  % 0‑indexed
    % D = min(max(D, 0), N-1);  % clip to equalizer span
    D = floor((N-1)/2);
end
