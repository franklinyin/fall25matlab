function D = pick_decision_delay(h, N)
% Choose a reasonable decision delay 0..N-1.
    % Older solution: place center tap near the largest‑energy channel tap.
    % [~, idx] = max(abs(h));
    % D = idx - 1;  % 0‑indexed
    % D = min(max(D, 0), N-1);  % clip to equalizer span

    % NEW solution:put the decision roughly in the middle of the equalizer
    D = floor((N-1)/2);
end
