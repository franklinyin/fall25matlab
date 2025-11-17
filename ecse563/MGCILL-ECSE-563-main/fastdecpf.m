function [V, delta, Psl, Qgv, N, time,Pf_MW, Qf_Mvar, Sf_MVA] = ...
    fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter,nfrom,nto)


tic;

% Setup initial parameters
nbus = length(Y);
G = real(Y);
B = imag(Y);

Vm = ones(nbus, 1);
Vm(ipv) = V0(ipv);
Vm(is)  = V0(is);
delta = zeros(nbus, 1);
V = Vm .* exp(1j * delta);

% Convert to per-unit
Pg = Pg / Sbase;
Qg = Qg / Sbase;
Pd = Pd / Sbase;
Qd = Qd / Sbase;

Pg_inj = Pg - Pd;
Qg_inj = Qg - Qd;

% Define index sets for pvpq
pvpq = setdiff([ipv(:); ipq(:)], is);  
pq   = setdiff(ipq(:), is);

% --- Build and scale B' and B'' matrices ---
Bp = -imag(Y);
Bq = -imag(Y);

% Remove slack
Bp = Bp(pvpq, pvpq);
Bq = Bq(pq, pq);

% Apply voltage scaling
Bp = Bp ./ (Vm(pvpq) * Vm(pvpq).');
Bq = Bq ./ (Vm(pq) * Vm(pq).');

% --- Iteration ---
for N = 1:maxiter
    Vm = abs(V);

    % Compute injected powers
    Pcalc = sum((Vm * Vm.') .* (G .* cos(delta - delta.') + B .* sin(delta - delta.')), 2);
    Qcalc = sum((Vm * Vm.') .* (G .* sin(delta - delta.') - B .* cos(delta - delta.')), 2);

    % Power mismatches
    dP = Pg_inj(pvpq) - Pcalc(pvpq);
    dQ = Qg_inj(pq)   - Qcalc(pq);

    % Convergence check
    if max(abs([dP; dQ])) < toler
        break;
    end

    % Fast-decoupled updates
    dDelta = Bp \ dP;
    dV     = Bq \ dQ;

    % Update
    delta(pvpq) = delta(pvpq) + dDelta;
    Vm(pq)      = Vm(pq) + dV;

    % Enforce fixed magnitudes
    Vm(ipv) = V0(ipv);
    Vm(is)  = V0(is);

    V = Vm .* exp(1j * delta);
end

time = toc;

Vcomplex = Vm .* exp(1j * delta);
I = Y * Vcomplex;

% Slack active power (MW)
Psl = real(Vcomplex(is) * conj(I(is))) * Sbase;

% Reactive power at all buses (Mvar)
Qgv = imag(Vcomplex .* conj(I)) * Sbase;

V = Vm;

% Get power flow per line
nl = length(nfrom);
Pf_MW  = zeros(nl,1);
Qf_Mvar= zeros(nl,1);
Sf_MVA = zeros(nl,1);

for k = 1:nl
    i = nfrom(k); 
    j = nto(k);
    Vi = Vcomplex(i); Vj = Vcomplex(j);

    y_series = -Y(i,j);
    Iij      = (Vi - Vj) * y_series;   
    Sij_pu   = Vi * conj(Iij);

    Pf_MW(k)   = real(Sij_pu) * Sbase;
    Qf_Mvar(k) = imag(Sij_pu) * Sbase;
    Sf_MVA(k)  = abs(Sij_pu)  * Sbase;


end

end
