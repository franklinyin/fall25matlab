function [h, descr] = generate_channel(ch)
% Returns an example ISI channel impulse response.
%   ch.type: 'low' | 'mild' | 'severe'

    if ~isfield(ch,'type'), ch.type='mild'; end
    switch lower(ch.type)
        case 'low'
            % Low-ISI real-valued 3-tap channel (no phase rotation)
            h = [1.0, 0.2, 0.05];
            descr = 'low-ISI 3-tap (real-valued)';
        case 'mild'
            % modest postcursor ISI
            h = [0.9 + 0.0j, 0.3 - 0.2j, 0.15 + 0.05j];
            descr = 'mild 3‑tap';
        case 'severe'
            % stronger multipath with precursor and postcursor
            h = [0.3 + 0.2j, 0.9 + 0.0j, 0.25 - 0.15j, 0.15 + 0.10j, 0.05 - 0.06j];
            descr = 'severe 5‑tap';
    end
end
