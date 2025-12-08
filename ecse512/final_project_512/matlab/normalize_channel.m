function h = normalize_channel(h)
% Normalize channel to unit energy
    h = h / norm(h);
end
