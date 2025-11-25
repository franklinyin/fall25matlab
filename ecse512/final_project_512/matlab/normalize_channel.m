function h = normalize_channel(h)
%NORMALIZE_CHANNEL Scales h[n] to unit energy.
    h = h / norm(h);
end
