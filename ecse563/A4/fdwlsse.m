
function [delta, V, Nit, tsec] = fdwlsse(nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)
%FDWLSSE Fast-decoupled weighted least-squares state estimator.
%   States: voltage angles at non-slack buses (delta) and |V| at PQ buses.
%   Measurements: P/Q injections, P/Q flows, V magnitudes with std devs.
%   All quantities are in per unit; line data include total line charging b.
%
%   Reference: ECSE 563 notes (Power Flow + State Estimation).

tic;
% Build Ybus
nb = max([nfrom(:); nto(:)]);
Y = zeros(nb, nb);
for k = 1:numel(nfrom)
    i = nfrom(k); j = nto(k);
    z = r(k) + 1j*x(k);
    y = 1/z;
    bc = 1j*b(k)/2;
    Y(i,i) = Y(i,i) + y + bc;
    Y(j,j) = Y(j,j) + y + bc;
    Y(i,j) = Y(i,j) - y;
    Y(j,i) = Y(j,i) - y;
end
G = real(Y); B = imag(Y);

% Extract measurement vectors and weights
zP = [Pinj(:,2); Pflow(:,3)];     sP = [Pinj(:,3); Pflow(:,4)];
zQ = [Qinj(:,2); Qflow(:,3)];     sQ = [Qinj(:,3); Qflow(:,4)];
zV = Vnode(:,2);                  sV = Vnode(:,3);
WP = diag(1./(sP.^2)); WQ = diag(1./(sQ.^2)); WV = diag(1./(sV.^2));

% Indexing
nbus = nb;
ref = 1;          % Slack bus index (assumed 1)
pv = [];          % Not used in FD decoupled SE (hold voltages via Vnode if any)
pq = setdiff(1:nbus, [ref, pv]);

% Initial state
delta = zeros(nbus,1); V = ones(nbus,1);
% If any Vnode provided, use them as starting magnitudes at those buses
if ~isempty(Vnode)
    V(Vnode(:,1)) = Vnode(:,2);
end

for Nit = 1:maxiter
    % --- P-subproblem: estimate delta (angles) ---
    % Calculate hP(delta) for injections and flows (approximation: use B matrix)
    % Injections: P_i ≈ sum_j -B_ij (delta_i - delta_j) for |V|≈1
    % Flows i->j : P_ij ≈ (delta_i - delta_j)/x_ij
    mPinj = size(Pinj,1); mPflow = size(Pflow,1);
    hPinj = zeros(mPinj,1); HP = zeros(mPinj + mPflow, nbus);
    % Injection rows
    for m = 1:mPinj
        i = Pinj(m,1);
        for j = 1:nbus
            if j==i, continue; end
            hPinj(m) = hPinj(m) + (-B(i,j))*(delta(i) - delta(j));
            HP(m,i) = HP(m,i) + (-B(i,j));
            HP(m,j) = HP(m,j) + (+B(i,j));
        end
    end
    % Flow rows
    hPflow = zeros(mPflow,1);
    for m = 1:mPflow
        i = Pflow(m,1); j = Pflow(m,2);
        % Find series reactance x_ij (assume unique line pair)
        idx = find((nfrom==i & nto==j) | (nfrom==j & nto==i), 1);
        xij = x(idx);
        hPflow(m) = (delta(i) - delta(j))/xij;
        HP(mPinj+m, i) = HP(mPinj+m, i) + 1/xij;
        HP(mPinj+m, j) = HP(mPinj+m, j) - 1/xij;
    end
    hP = [hPinj; hPflow];
    % Remove reference angle column/row (delta_ref = 0)
    keep = setdiff(1:nbus, ref);
    HPk = HP(:, keep);
    rP = [zP] - hP;
    % Solve normal equations
    d_delta = (HPk.' * WP * HPk) \ (HPk.' * WP * rP);
    delta(keep) = delta(keep) + d_delta;

    % --- Q/V-subproblem: estimate |V| (magnitudes) ---
    % Measurements: Q injections, Q flows, |V| magnitudes
    mQinj = size(Qinj,1); mQflow = size(Qflow,1); mV = size(Vnode,1);
    hQinj = zeros(mQinj,1); HQ = zeros(mQinj + mQflow + mV, numel(pq));
    % Approximation: use susceptance-looking Jacobian (fast-decoupled)
    for m = 1:mQinj
        i = Qinj(m,1);
        % Q_i ≈ -B_ii*(V_i - 1) - sum_{j≠i} B_ij*(V_j - 1)  (linearized)
        hQinj(m) = -B(i,i)*(V(i)-1);
        for j = 1:nbus
            if j==i, continue; end
            hQinj(m) = hQinj(m) - B(i,j)*(V(j)-1);
        end
        % Jacobian wrt V magnitudes at PQ buses
        for k = 1:numel(pq)
            j = pq(k);
            if j==i
                HQ(m,k) = HQ(m,k) - B(i,i);
            else
                HQ(m,k) = HQ(m,k) - B(i,j);
            end
        end
    end
    % Q flows and V magnitudes
    hQflow = zeros(mQflow,1); HV = zeros(mV, numel(pq));
    for m = 1:mQflow
        i = Qflow(m,1); j = Qflow(m,2);
        idx = find((nfrom==i & nto==j) | (nfrom==j & nto==i), 1);
        bij = -1/x(idx);  %#ok<NASGU>  % not used in very rough decoupled Q model
        hQflow(m) = 0;    % neglected in the fast-decoupled linearization
    end
    hV = V(Vnode(:,1));
    for m = 1:mV
        k = find(pq==Vnode(m,1));
        if ~isempty(k), HV(m,k) = 1; end
    end
    HQbig = [HQ; zeros(mQflow, size(HQ,2)); HV];
    hQ = [hQinj; hQflow; hV];

    rQ = [zQ; zeros(mQflow,1); zV] - hQ;
    dV = (HQbig.' * blkdiag(WQ, eye(mQflow), WV) * HQbig) \ ...
         (HQbig.' * blkdiag(WQ, eye(mQflow), WV) * rQ);
    V(pq) = V(pq) + dV;

    % Convergence (use infinity norm of updates)
    if max(abs(d_delta)) < toler && (isempty(dV) || max(abs(dV)) < toler)
        break;
    end
end

tsec = toc;
end
