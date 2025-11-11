clear; clc;

% (Eqs. 7.28a–b)
Ap_amp = 0.89125;
As_amp = 0.17783;
wp = 0.2*pi;
ws = 0.3*pi;
Td = 1;
Rp = -20*log10(Ap_amp);
Rs = -20*log10(As_amp);

% Prewarp the edges (Eq. 7.26 → 7.29)
Wp = (2/Td)*tan(wp/2);
Ws = (2/Td)*tan(ws/2);

% Order from Eq. (7.33); cutoff from (7.32b)

N_real = log( ((1/As_amp)^2 - 1) / ((1/Ap_amp)^2 - 1) ) ...
         / ( 2*log(Ws/Wp) );
N = ceil(N_real);

% choose Omega_c to meet stopband exactly (book does this): Eq. (7.32b)
Oc = Ws / (((1/As_amp)^2 - 1)^(1/(2*N)));  % analog 3 dB cutoff

fprintf('Eq. (7.33): N* = %.6f  → choose N = %d\n', N_real, N);

%continuous-time Butterworth Hc(s)
k  = 1:(N/2);
a1 = 2*Oc*sin((2*k-1)*pi/(2*N));
a0 = Oc^2 * ones(size(k));

Kc = Oc^N; % (so that Hc(0)=1)

% build full denominator polynomial from biquads
% Each biquad is (s^2 + a1*s + a0)
as = [1];  % start with 1
for i = 1:numel(k)
    biquad = [1, a1(i), a0(i)];
    as = conv(as, biquad);
end

if mod(N, 2) == 1
    as = conv(as, [1, Oc]);
end

bs = Kc;


% bilinear transform to get H(z) (Td=1)
Td = 1;
[bz, az] = manual_bilinear(bs, as, Td);

% Normalize to ensure DC gain is 1 (H(1) = 1)
% At DC, z=1, so H(1) = sum(bz) / sum(az)
dc_gain = sum(bz) / sum(az);
if abs(dc_gain - 1.0) > 1e-10
    bz = bz / dc_gain;  % normalize numerator to get H(1) = 1
    fprintf('Normalized filter: DC gain was %.6f, normalized to 1.0\n', dc_gain);
end


A = zeros(numel(k), 3);
for i = 1:numel(k)
    den_i = 4 + 2*a1(i) + a0(i); % normalizing constant
    A(i,:) = [1, (-8 + 2*a0(i))/den_i, (4 - 2*a1(i) + a0(i))/den_i];
end

g = Kc / prod(4 + 2*a1 + a0);
p = 1;  % build (1 + z^-1)^N
bin = 1;
for i = 1:N
    bin = conv(bin, [1 1]);
end
b_factored = g * bin;

% print (7.35)-style result
fprintf('H(z) (Eq. 7.35 form):\n');
fprintf('  Numerator: %.7f * (1 + z^{-1})^%d\n', g, N);
for i = 1:size(A,1)
    fprintf('  Denominator section %d: 1 %+.4f z^{-1} %+0.4f z^{-2}\n', ...
        i, A(i,2), A(i,3));
end

%---------------------------------------
% Fig. 7.11-style plots
%---------------------------------------
% Check DC gain (should be 1)
z_dc = 1;  % at ω=0, z = exp(j*0) = 1
H_dc = manual_freqz(bz, az, z_dc);
fprintf('DC gain check: H(1) = %.6f (should be 1.0)\n', H_dc);
fprintf('  sum(bz) = %.6f, sum(az) = %.6f\n', sum(bz), sum(az));
fprintf('  sum(bz)/sum(az) = %.6f\n\n', sum(bz)/sum(az));

nfft = 4096;
% Manual frequency response calculation
w = linspace(0, pi, nfft)';
z = exp(1j*w);  % z = e^(jω)
H = manual_freqz(bz, az, z);
mag = abs(H);
magdB = 20*log10(mag);

% Group delay: -d(angle(H))/dω
wgd = w;
gd = manual_grpdelay(bz, az, z, w);

% Figure 1: Log magnitude
figure('Name','Log magnitude','Color','w');
plot(w/pi, magdB, 'LineWidth', 1.2); grid on;
xline(0.2,'--'); xline(0.3,'--');
ylabel('Magnitude (dB)'); xlabel('\omega/\pi'); 
xlim([0 1]); ylim([-100 20]);
xticks([0.2 0.4 0.6 0.8 1.0]);
xticklabels({'0.2\pi', '0.4\pi', '0.6\pi', '0.8\pi', '\pi'});
title('Fig. 7.11(a): Log magnitude');
saveas(gcf, 'log_magnitude.png');

% Figure 2: Magnitude
figure('Name','Magnitude','Color','w');
plot(w/pi, mag, 'LineWidth', 1.2); grid on;
xline(0.2,'--'); xline(0.3,'--');
ylabel('Magnitude'); xlabel('\omega/\pi'); 
xlim([0 1]); ylim([0 1.2]);
xticks([0.2 0.4 0.6 0.8 1.0]);
xticklabels({'0.2\pi', '0.4\pi', '0.6\pi', '0.8\pi', '\pi'});
title('Fig. 7.11(b): Magnitude');
saveas(gcf, 'magnitude.png');

% Figure 3: Group delay
figure('Name','Group delay','Color','w');
plot(wgd/pi, gd, 'LineWidth', 1.2); grid on;
xline(0.2,'--'); xline(0.3,'--');
xlabel('\omega/\pi'); ylabel('Group delay [samples]'); 
xlim([0 1]); ylim([0 12]);
xticks([0.2 0.4 0.6 0.8 1.0]);
xticklabels({'0.2\pi', '0.4\pi', '0.6\pi', '0.8\pi', '\pi'});
title('Fig. 7.11(c): Group delay');
saveas(gcf, 'group_delay.png');



%% helper functions

function [bz, az] = manual_bilinear(bs, as, Td)
    if as(1) ~= 1
        bs = bs / as(1);
        as = as / as(1);
    end
    
    M = length(bs);
    N = length(as);
    
    T = 2/Td;
    
    one_plus_z = [1, 1];
    one_plus_z_power = 1;
    for i = 1:(N-1)
        one_plus_z_power = conv(one_plus_z_power, one_plus_z);
    end
    
    bz_temp = zeros(1, length(one_plus_z_power));
    for k = 1:M
        if bs(k) ~= 0
            power_s = k - 1; 
            one_minus_z = [1, -1];  % (1 - z^-1)
            one_minus_z_power = 1;
            for i = 1:power_s
                one_minus_z_power = conv(one_minus_z_power, one_minus_z);
            end
            
            one_plus_z_power2 = 1;
            if (N-1-power_s) > 0
                for i = 1:(N-1-power_s)
                    one_plus_z_power2 = conv(one_plus_z_power2, one_plus_z);
                end
            end
            
            term = bs(k) * (T^power_s) * conv(one_minus_z_power, one_plus_z_power2);
            
            if length(term) < length(bz_temp)
                term = [term, zeros(1, length(bz_temp) - length(term))];
            elseif length(term) > length(bz_temp)
                bz_temp = [bz_temp, zeros(1, length(term) - length(bz_temp))];
            end
            bz_temp = bz_temp + term;
        end
    end
    
    az_temp = zeros(1, length(one_plus_z_power));
    for k = 1:N
        if as(k) ~= 0
            power_s = k - 1;
            one_minus_z = [1, -1]; 
            one_minus_z_power = 1;
            for i = 1:power_s
                one_minus_z_power = conv(one_minus_z_power, one_minus_z);
            end
            
            one_plus_z_power2 = 1;
            if (N-1-power_s) > 0
                for i = 1:(N-1-power_s)
                    one_plus_z_power2 = conv(one_plus_z_power2, one_plus_z);
                end
            end
            
            term = as(k) * (T^power_s) * conv(one_minus_z_power, one_plus_z_power2);
            
            % pad to match length
            if length(term) < length(az_temp)
                term = [term, zeros(1, length(az_temp) - length(term))];
            elseif length(term) > length(az_temp)
                az_temp = [az_temp, zeros(1, length(term) - length(az_temp))];
            end
            az_temp = az_temp + term;
        end
    end
    
    bz = bz_temp / az_temp(1);
    az = az_temp / az_temp(1);
end

function H = manual_freqz(bz, az, z)
    % Manual frequency response: H(z) = B(z)/A(z) evaluated at z values
    % bz, az are polynomial coefficients in z^-1 (standard MATLAB format)
    % z can be scalar, vector, or array
    
    % Evaluate numerator B(z) = sum(bz(k) * z^-(k-1))
    % Use explicit computation: z^-n = 1/(z^n) for better numerical stability
    B = zeros(size(z));
    for k = 1:length(bz)
        if bz(k) ~= 0
            n = k - 1;  % power of z^-1
            if n == 0
                B = B + bz(k);
            else
                B = B + bz(k) ./ (z.^n);
            end
        end
    end
    
    % Evaluate denominator A(z) = sum(az(k) * z^-(k-1))
    A = zeros(size(z));
    for k = 1:length(az)
        if az(k) ~= 0
            n = k - 1;  % power of z^-1
            if n == 0
                A = A + az(k);
            else
                A = A + az(k) ./ (z.^n);
            end
        end
    end
    
    H = B ./ A;
end

function gd = manual_grpdelay(bz, az, z, w)
    
    H = manual_freqz(bz, az, z);
    phase = angle(H);
    
    gd = zeros(size(w));
    if length(w) > 1
        gd(2:end-1) = -(phase(3:end) - phase(1:end-2)) ./ (w(3:end) - w(1:end-2));
        gd(1) = -(phase(2) - phase(1)) / (w(2) - w(1));
        gd(end) = -(phase(end) - phase(end-1)) / (w(end) - w(end-1));
    end
end


