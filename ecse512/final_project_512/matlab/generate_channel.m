function [h, descr] = generate_channel(ch)
%GENERATE_CHANNEL Returns an example ISI channel impulse response.
%   ch.type: 'mild' | 'severe' | 'random'
%   ch.K   : length when 'random'
%   ch.randSTD : std for Rayleigh taps (before normalization)

    if ~isfield(ch,'type'), ch.type='mild'; end
    switch lower(ch.type)
        case 'mild'
            % modest postcursor ISI
            h = [0.9 + 0.0j, 0.3 - 0.2j, 0.15 + 0.05j];
            descr = 'mild 3‑tap';
        case 'severe'
            % stronger multipath with precursor and postcursor
            h = [0.3 + 0.2j, 0.9 + 0.0j, 0.25 - 0.15j, 0.15 + 0.10j, 0.05 - 0.06j];
            descr = 'severe 5‑tap';
        otherwise  % 'random'
            if ~isfield(ch,'K'), ch.K = 5; end
            if ~isfield(ch,'randSTD'), ch.randSTD = 1; end
            h = ch.randSTD/sqrt(2) * (randn(1,ch.K) + 1j*randn(1,ch.K));
            descr = sprintf('random Rayleigh %d‑tap', ch.K);
    end
end
