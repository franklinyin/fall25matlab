% Q4 implementation (fast-decoupled WLS SE)
function [delta, V, Niter, elapsed] = fdwlsse(nfrom, nto, r, x, b, ...
    Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)

    % Number of buses
    nbus = max(max(nfrom), max(nto));

    % Network admittance matrix
    Ybus = admittance(nfrom, nto, r, x, b);

    % Measurements and weights
    [z, W] = extract_measurements(Pinj, Qinj, Pflow, Qflow, Vnode);

    % Flat start state: [delta(2..N); V(1..N)]
    x_state = [zeros(nbus-1,1); ones(nbus,1)];

    % ---------------------------------------------------------------------
    % Fast‑decoupled WLS setup:
    %   - Partition measurements into active (P, Pf), reactive (Q, Qf),
    %     and voltage magnitude (V)
    %   - Build decoupled Jacobian blocks at flat profile
    %   - Precompute gain matrices G_AA (angle block) and G_RR (voltage block)
    % ---------------------------------------------------------------------
    [idxA, idxR, idxV] = build_measurement_indices(Pinj, Qinj, Pflow, Qflow, Vnode);
    idxA  = idxA(:);            % active P/Pf measurements
    idxR  = idxR(:);            % reactive Q/Qf measurements
    idxV  = idxV(:);            % V measurements
    idxQV = [idxR; idxV];       % reactive + voltage block

    % Jacobian at flat profile (used for the whole iteration = "fast" part)
    h0  = compute_measurements(x_state, Ybus, nbus, nfrom, nto, r, x, b, ...
                               Pinj, Qinj, Pflow, Qflow, Vnode);
    H0  = compute_jacobian(x_state, h0, Ybus, nbus, nfrom, nto, r, x, b, ...
                           Pinj, Qinj, Pflow, Qflow, Vnode);

    % Split into angle and magnitude sensitivities
    Hdelta = H0(:, 1:nbus-1);   % ∂h/∂δ
    HV     = H0(:, nbus:end);   % ∂h/∂V

    % Active measurements use only ∂P/∂δ, ∂Pf/∂δ (P−δ block)
    H_A = Hdelta(idxA, :);
    W_A = W(idxA, idxA);
    G_AA = H_A.' * (W_A * H_A);

    % Reactive + voltage measurements use only ∂Q/∂V, ∂Qf/∂V and dV/dV (Q−V block)
    H_R = HV(idxQV, :);
    W_R = W(idxQV, idxQV);
    G_RR = H_R.' * (W_R * H_R);

    % ---------------------------------------------------------------------
    % Fast‑decoupled WLS iterations
    % ---------------------------------------------------------------------
    tic;
    for k = 1:maxiter
        % Nonlinear measurement prediction at current estimate
        h = compute_measurements(x_state, Ybus, nbus, nfrom, nto, r, x, b, ...
                                 Pinj, Qinj, Pflow, Qflow, Vnode);
        residual = z - h;   % r = z − h(x^k)

        % Partition residuals
        r_A  = residual(idxA);      % active P/Pf residuals
        r_RV = residual(idxQV);     % reactive Q/Qf + voltage residuals

        % Decoupled right‑hand sides
        t_A = H_A.' * (W_A * r_A);
        t_R = H_R.' * (W_R * r_RV);

        % Solve the two smaller linear systems:
        %   G_AA * Δδ = t_A
        %   G_RR * ΔV = t_R
        d_delta = G_AA \ t_A;
        d_V     = G_RR \ t_R;

        % Assemble full state increment and update
        dx = [d_delta; d_V];
        x_state = x_state + dx;

        % Convergence test on state change (same criterion as before)
        if max(abs(dx)) < toler
            break;
        end
    end
    elapsed = toc;
    Niter   = k;

    % Recover full angle vector (slack = 0) and voltage magnitudes
    delta = [0; x_state(1:nbus-1)];
    V     = x_state(nbus:end);
end

% -------------------------------------------------------------------------
% For the convenience of immediate reference, the helper functions are
% included underneath the body of the code
% -------------------------------------------------------------------------

%% Helper: build measurement index sets for decoupling
function [idxA, idxR, idxV] = build_measurement_indices(Pinj, Qinj, Pflow, Qflow, Vnode)
    k = 0;
    idxA = [];
    idxR = [];
    idxV = [];

    if ~isempty(Pinj)
        nP = size(Pinj,1);
        idxA = [idxA; (k+1 : k+nP).'];
        k = k + nP;
    end
    if ~isempty(Qinj)
        nQ = size(Qinj,1);
        idxR = [idxR; (k+1 : k+nQ).'];
        k = k + nQ;
    end
    if ~isempty(Pflow)
        nPf = size(Pflow,1);
        idxA = [idxA; (k+1 : k+nPf).'];
        k = k + nPf;
    end
    if ~isempty(Qflow)
        nQf = size(Qflow,1);
        idxR = [idxR; (k+1 : k+nQf).'];
        k = k + nQf;
    end
    if ~isempty(Vnode)
        nV = size(Vnode,1);
        idxV = (k+1 : k+nV).';
        k = k + nV;
    end
end

%% Extract measurements and build weight matrix
function [z, W] = extract_measurements(Pinj, Qinj, Pflow, Qflow, Vnode)
    z = [];
    sigma2 = [];
    
    if ~isempty(Pinj)
        z = [z; Pinj(:,2)];
        sigma2 = [sigma2; Pinj(:,3).^2];
    end
    if ~isempty(Qinj)
        z = [z; Qinj(:,2)];
        sigma2 = [sigma2; Qinj(:,3).^2];
    end
    if ~isempty(Pflow)
        z = [z; Pflow(:,3)];
        sigma2 = [sigma2; Pflow(:,4).^2];
    end
    if ~isempty(Qflow)
        z = [z; Qflow(:,3)];
        sigma2 = [sigma2; Qflow(:,4).^2];
    end
    if ~isempty(Vnode)
        z = [z; Vnode(:,2)];
        sigma2 = [sigma2; Vnode(:,3).^2];
    end
    
    W = diag(1 ./ sigma2);
end

%% Compute predicted measurements
function h = compute_measurements(x, Ybus, nbus, nfrom, nto, r, x_param, b, ...
                                  Pinj, Qinj, Pflow, Qflow, Vnode)
    
    delta = [0; x(1:nbus-1)];
    V = x(nbus:end);
    Vcmp = V .* exp(1j*delta);
    
    [Pbus, Qbus] = compute_injections(Vcmp, Ybus);
    
    h = [];
    if ~isempty(Pinj)
        h = [h; Pbus(Pinj(:,1))];
    end
    if ~isempty(Qinj)
        h = [h; Qbus(Qinj(:,1))];
    end
    if ~isempty(Pflow)
        h = [h; compute_pflow(Pflow, Vcmp, nfrom, nto, r, x_param, b)];
    end
    if ~isempty(Qflow)
        h = [h; compute_qflow(Qflow, Vcmp, nfrom, nto, r, x_param, b)];
    end
    if ~isempty(Vnode)
        h = [h; abs(Vcmp(Vnode(:,1)))];
    end
end

%% Compute bus injections
function [Pbus, Qbus] = compute_injections(Vcmp, Ybus)
    Iinj = Ybus * Vcmp;
    Sbus = Vcmp .* conj(Iinj);
    Pbus = real(Sbus);
    Qbus = imag(Sbus);
end

%% Compute active power flows
function Pflow_vals = compute_pflow(Pflow, Vcmp, nfrom, nto, r, x, b)
    Pflow_vals = zeros(size(Pflow,1),1);
    yser = 1 ./ (r + 1j*x);
    bsh = 1j * b / 2;
    
    for k = 1:size(Pflow,1)
        i = Pflow(k,1);
        j = Pflow(k,2);
        [idx, i2, j2] = find_line_index(i, j, nfrom, nto);
        y = yser(idx);
        bs = bsh(idx);
        Iij = (Vcmp(i2) - Vcmp(j2)) * y + Vcmp(i2) * bs;
        Sij = Vcmp(i2) * conj(Iij);
        Pflow_vals(k) = real(Sij);
    end
end

%% Compute reactive power flows
function Qflow_vals = compute_qflow(Qflow, Vcmp, nfrom, nto, r, x, b)
    Qflow_vals = zeros(size(Qflow,1),1);
    yser = 1 ./ (r + 1j*x);
    bsh = 1j * b / 2;
    
    for k = 1:size(Qflow,1)
        i = Qflow(k,1);
        j = Qflow(k,2);
        [idx, i2, j2] = find_line_index(i, j, nfrom, nto);
        y = yser(idx);
        bs = bsh(idx);
        Iij = (Vcmp(i2) - Vcmp(j2)) * y + Vcmp(i2) * bs;
        Sij = Vcmp(i2) * conj(Iij);
        Qflow_vals(k) = imag(Sij);
    end
end

%% Find line index and direction
function [idx, i2, j2] = find_line_index(i, j, nfrom, nto)
    idx = find(nfrom==i & nto==j, 1);
    if isempty(idx)
        idx = find(nfrom==j & nto==i, 1);
        i2 = j;
        j2 = i;
    else
        i2 = i;
        j2 = j;
    end
end

%% Compute numerical Jacobian
function H = compute_jacobian(x, h0, Ybus, nbus, nfrom, nto, r, x_param, b, ...
                               Pinj, Qinj, Pflow, Qflow, Vnode)
    
    nstate = numel(x);
    m = numel(h0);
    H = zeros(m, nstate);
    eps_fd = 1e-6;
    
    for j = 1:nstate
        xpert = x;
        xpert(j) = xpert(j) + eps_fd;
        hpert = compute_measurements(xpert, Ybus, nbus, nfrom, nto, r, x_param, b, ...
                                     Pinj, Qinj, Pflow, Qflow, Vnode);
        H(:,j) = (hpert - h0) / eps_fd;
    end
end
