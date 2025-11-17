function [V, delta, Psl, Qgv, N, time,Pf_MW,Qf_Mvar,Sf_MVA] = ... 
    nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter, nfrom, nto)

% Start timer and track iterations
tic
N = 0;
nbus = length(Y);

% Initialize V maginitude and delta, and setup ipv  
Vm = ones(nbus,1);
Vm(ipv) = V0(ipv);
delta = zeros(nbus,1);
V = Vm .* exp(1j*delta);

% Convert power to common base
Pg = Pg/Sbase;
Qg = Qg/Sbase;   
Pd = Pd/Sbase;   
Qd = Qd/Sbase;   

% Find the net power per node (generated - demanded)
Pg_inj = Pg - Pd;
Qg_inj = Qg - Qd;

% Y = G + jB
G = real(Y); 
B = imag(Y);

pvpq = setdiff([ipv(:); ipq(:)], is);
pq   = setdiff(ipq(:), is);


for N = 1:maxiter

    Vm = abs(V);

    Pcalc = sum( (Vm * Vm.') .* ( G .* cos(delta - delta.') + B .* sin(delta - delta.') ), 2 );
    Qcalc = sum( (Vm * Vm.') .* ( G .* sin(delta - delta.') - B .* cos(delta - delta.') ), 2 );

    % find the difference, seperated feed give P and Q into into Jacobian
    delta_P = Pg_inj - Pcalc;
    delta_Q = Qg_inj - Qcalc;


    mismatch = [ delta_P(pvpq); delta_Q(pq) ];


    % Convergence test
    if max(abs(mismatch)) < toler
        break
    end

    % Get jacobian
    [J, H_jac, N_jac, M_jac, L_jac] = buildJacobian(G, B, Vm, delta, Pcalc, Qcalc);

    % Reduce jacobian to include PV and PQ indices of intrest
    Jred = [ H_jac(pvpq,pvpq),  N_jac(pvpq, pq);
         M_jac(pq, pvpq), L_jac(pq, pq) ];


    dx = Jred \ mismatch;

    d_delta = dx(1:numel(pvpq));
    d_V     = dx(numel(pvpq)+1:end) .*  abs(V(pq));

    
    % Update state
    delta(pvpq) = delta(pvpq) + d_delta;  % update angles for PV+PQ
    Vm(pq)      = Vm(pq)      + d_V;      % update magnitudes for PQ only

    % % 
    % Vm(ipv) = V0(ipv);
    % Vm(is)  = V0(is);

    
    V = Vm .* exp(1j*delta);

   
end


Vcomplex = Vm .* exp(1j*delta);
I = Y * Vcomplex;

% Slack active power (MW)
Psl = real( Vcomplex(is) * conj(I(is)) ) * Sbase;

% Reactive power at PV buses (Mvar) in same order as ipv
Qgv = imag(Vcomplex .* conj(I)) * Sbase;

V     = Vm;
time  = toc;

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


function [J, H, N, M, L] = buildJacobian(G, B, V, delta, P, Q)

nb = length(V);


% --- allocate blocks (all buses, full size) ---
H = zeros(nb, nb);   % ∂P/∂δ
N = zeros(nb, nb);   % ∂P/∂|V|
M = zeros(nb, nb);   % ∂Q/∂δ
L = zeros(nb, nb);   % ∂Q/∂|V|

% --- fill entries using the table ---
for i = 1:nb
    for j = 1:nb
        if i == j
            H(i,j) = -Q(i) - B(i,i)*V(i)^2;
            N(i,j) =  P(i) + G(i,i)*V(i)^2;
            M(i,j) =  P(i) - G(i,i)*V(i)^2;
            L(i,j) =  Q(i) - B(i,i)*V(i)^2;
        else
            dij = delta(i) - delta(j);
            Hv = V(i)*V(j)*( G(i,j)*sin(dij) - B(i,j)*cos(dij) );
            Nv = V(i)*V(j)*( G(i,j)*cos(dij) + B(i,j)*sin(dij) );
            H(i,j) = Hv;
            N(i,j) = Nv;
            M(i,j) = -Nv;
            L(i,j) = Hv;
        end
    end
end

J = [H N; M L];
end
