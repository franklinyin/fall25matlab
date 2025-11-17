function [V, delta, Psl, Qgv, N, time,Pf_MW, Qf_Mvar, Sf_MVA] = ... 
    decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter,nfrom,nto)

tic;

nbus = length(Y);
G = real(Y);
B = imag(Y);

% Setup initial parameters
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

% Build constant Jacobian blocks
Bp = -imag(Y);
Bq = -imag(Y);

% remove slack row/column
Bp = Bp(pvpq, pvpq);
Bq = Bq(pq, pq);


for N = 1:maxiter
    Vm = abs(V);

    % Compute injected powers
    Pcalc = sum( (Vm * Vm.') .* (G .* cos(delta - delta.') + B .* sin(delta - delta.')), 2 );
    Qcalc = sum( (Vm * Vm.') .* (G .* sin(delta - delta.') - B .* cos(delta - delta.')), 2 );

    % get power mismatches
    dP = Pg_inj(pvpq) - Pcalc(pvpq);
    dQ = Qg_inj(pq)   - Qcalc(pq);

    % Convergence check
    if max(abs([dP; dQ])) < toler
        break;
    end

    % Decoupled updates (normalized by voltage magnitudes)
    dDelta = Bp \ (dP ./ Vm(pvpq));
    dV     = Bq \ (dQ ./ Vm(pq));

    % Update state variables
    delta(pvpq) = delta(pvpq) + dDelta;
    Vm(pq)      = Vm(pq) + dV;

    % Enforce PV and Slack voltage magnitudes
    Vm(ipv) = V0(ipv);
    Vm(is)  = V0(is);

    % Recompute complex voltage
    V = Vm .* exp(1j * delta);
end

time = toc;

Vcomplex = Vm .* exp(1j * delta);
I = Y * Vcomplex;

% Slack active power
Psl = real(Vcomplex(is) * conj(I(is))) * Sbase;

% Reactive power injections at all buses 
Qgv = imag(Vcomplex .* conj(I)) * Sbase;

V = Vm;   % return magnitudes

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
