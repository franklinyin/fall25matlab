function h = normalize_channel(h)
% Scales h[n] to unit energy.
    h = h / norm(h);
end
