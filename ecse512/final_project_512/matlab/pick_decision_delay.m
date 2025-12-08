function D = pick_decision_delay(h, N)
% Choose decision delay - place near largest channel tap
    [~, idx] = max(abs(h));
    D = idx - 1;
    D = min(max(D, 0), N-1);
end
