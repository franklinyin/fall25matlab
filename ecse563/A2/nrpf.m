function [V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    % nrpf - Newton-Raphson AC power flow (polar coordinates)
    % Input:
    %   Y - admittance matrix
    %   is - slack bus index
    %   ipq, ipv - PQ and PV bus indices
    %   Pg, Qg - generation vectors
    %   Pd, Qd - demand vectors
    %   V0 - initial voltage magnitudes
    %   Sbase - base power (MVA)
    %   toler, maxiter - convergence parameters
    % Output:
    %   V - voltage magnitudes (p.u.)
    %   delta - voltage angles (radians)
    %   Psl - slack bus active power (MW)
    %   Qgv - PV bus reactive power (MVAr)
    %   N - number of iterations
    %   time - computation time (seconds)
%NRPF  Newton–Raphson AC power flow (polar coordinates).
%   [V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%
% Inputs:
%   Y      : (n x n) complex bus-admittance matrix
%   is     : slack bus index (1-based)
%   ipq    : vector of PQ bus indices (1-based)
%   ipv    : vector of PV bus indices (1-based); may include 'is' (removed internally)
%   Pg,Qg  : (n x 1) specified generation [MW, MVAr]
%   Pd,Qd  : (n x 1) specified demand    [MW, MVAr]
%   V0     : specified |V| (p.u.) for generator buses in 'ipv' order
%   Sbase  : MVA base
%   toler  : convergence tolerance on mismatch infinity-norm (p.u.)
%   maxiter: maximum allowed Newton iterations
%
% Outputs:
%   V      : (n x 1) bus voltage magnitudes [p.u.]
%   delta  : (n x 1) bus voltage angles [rad] (delta(is)=0 reference)
%   Psl    : slack active generation [MW]
%   Qgv    : (n x 1) reactive generation at PV buses [MVAr], zeros elsewhere
%   N      : iterations to convergence
%   time   : CPU time [s]
%
% Notes:
% - No matrix inversion is used; linear solves use backslash.
% - Jacobian blocks H,N,M,L follow standard NR polar formulation (course slides).
%   See ECSE563 Power Flow notes: NR derivation, mismatch definitions, and update rule. 
%   :contentReference[oaicite:10]{index=10} :contentReference[oaicite:11]{index=11} :contentReference[oaicite:12]{index=12}
%
n = size(Y,1);
G = real(Y); B = imag(Y);

% sets (ensure column vectors, unique, sorted)
ipq = unique(ipq(:));
ipv = setdiff(unique(ipv(:)), is);   % remove slack if present

% per-unit specified injections
Pinj = (Pg - Pd)/Sbase;   % p.u.
Qinj = (Qg - Qd)/Sbase;   % p.u.

% Initial guess
V     = ones(n,1);
delta = zeros(n,1);

% Apply specified |V| to generator buses (V0 is given in order of 'ipv' in data file).
% If user provided V0 including the slack in 'ipv', map by indices safely.
if numel(V0) == numel(ipv)
    V(ipv) = V0(:);
else
    % If V0 was provided for [is; ipv] as in the sample file (three generators),
    % map the entries by matching indices.
    Vgen_idx = [is; ipv(:)];
    Vgen_idx = unique(Vgen_idx(:),'stable');
    k = min(numel(V0), numel(Vgen_idx));
    V(Vgen_idx(1:k)) = V0(1:k);
end

% Indexing for unknowns (order: all non-slack angles [PV;PQ], then PQ magnitudes)
iang = [ipv; ipq];                    % angles unknown at all non-slack
iV   = ipq;                           % magnitudes unknown at PQ only

tic
for N = 1:maxiter
    % Calculate nodal P,Q from current (V,delta)
    P = zeros(n,1); Q = zeros(n,1);
    for i = 1:n
        for k = 1:n
            P(i) = P(i) + V(i)*V(k)*( G(i,k)*cos(delta(i)-delta(k)) + B(i,k)*sin(delta(i)-delta(k)) );
            Q(i) = Q(i) + V(i)*V(k)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
        end
    end

    % Mismatches (implicit equations only)
    dP = Pinj - P;             % all buses, but slack eq. not used
    dQ = Qinj - Q;             % all buses, only PQ used

    F  = [ dP(iang) ; dQ(iV) ];     % stacked mismatch (order must match Jacobian)

    % Convergence check
    if norm(F, inf) < toler
        break
    end

    % Build Jacobian sub-blocks H,N,M,L (see course notes) :contentReference[oaicite:13]{index=13}
    np = numel(iang); nq = numel(iV);
    H = zeros(np,np); Nblk = zeros(np,nq); M = zeros(nq,np); L = zeros(nq,nq);

    for a = 1:np
        i = iang(a);
        for b = 1:np
            k = iang(b);
            if i == k
                H(a,b) = -Q(i) - V(i)^2*B(i,i);
            else
                H(a,b) = V(i)*V(k)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
            end
        end
        for b = 1:nq
            k = iV(b);
            if i == k
                Nblk(a,b) = P(i)/V(i) + G(i,i)*V(i);
            else
                Nblk(a,b) = V(i)*( G(i,k)*cos(delta(i)-delta(k)) + B(i,k)*sin(delta(i)-delta(k)) );
            end
        end
    end

    for a = 1:nq
        i = iV(a);
        for b = 1:np
            k = iang(b);
            if i == k
                M(a,b) =  P(i) - V(i)^2*G(i,i);
            else
                M(a,b) = -V(i)*V(k)*( G(i,k)*cos(delta(i)-delta(k)) + B(i,k)*sin(delta(i)-delta(k)) );
            end
        end
        for b = 1:nq
            k = iV(b);
            if i == k
                L(a,b) = Q(i)/V(i) - B(i,i)*V(i);
            else
                L(a,b) = V(i)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
            end
        end
    end

    % Solve for updates and update state (no inv) :contentReference[oaicite:14]{index=14}
    J = [H Nblk; M L];
    dx = J \ F;

    ddelta = dx(1:np);
    dV     = dx(np+1:end);

    delta(iang) = delta(iang) + ddelta;
    V(iV)       = V(iV)       + dV;

    % Enforce slack angle reference and PV magnitudes
    delta(is) = 0;
    % Keep PV magnitudes at specified setpoints
    if numel(V0) == numel(ipv), V(ipv) = V0(:); end
end
time = toc;

% Final quantities
% Slack P generation (MW):  Psl = Pcalc_slack + Pd_slack - Pg_slack (all in p.u., then ×Sbase)
Psl = ( P(is) + (Pd(is)-Pg(is))/Sbase ) * Sbase;

% Reactive generation at PV buses (MVAr): Qg_i = Qcalc_i + Qd_i  (in p.u. × Sbase)
Qgv = zeros(n,1);
for k = 1:numel(ipv)
    i = ipv(k);
    Qgv(i) = ( Q(i) + Qd(i)/Sbase ) * Sbase;
end

% If we broke due to convergence early, N is the number of completed iterations.
if norm([ dP(iang) ; dQ(iV) ], inf) >= toler && N == maxiter
    warning('NRPF:MaxIter','NR did not meet tolerance in maxiter.');
end
end
