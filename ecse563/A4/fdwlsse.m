% Q4 implementation
function [delta, V, Niter, elapsed] = fdwlsse(nfrom, nto, r, x, b, ...
    Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)

    nbus = max(max(nfrom), max(nto));
    Ybus = admittance(nfrom, nto, r, x, b);
    [z, W] = extract_measurements(Pinj, Qinj, Pflow, Qflow, Vnode);
    
    x_state = [zeros(nbus-1,1); ones(nbus,1)];
    
    tic;
    for k = 1:maxiter
        h = compute_measurements(x_state, Ybus, nbus, nfrom, nto, r, x, b, ...
                                Pinj, Qinj, Pflow, Qflow, Vnode);
        residual = z - h;
        
        H = compute_jacobian(x_state, h, Ybus, nbus, nfrom, nto, r, x, b, ...
                            Pinj, Qinj, Pflow, Qflow, Vnode);
        
        G = H.' * (W * H);
        dx = G \ (H.' * (W * residual));
        x_state = x_state + dx;
        
        if max(abs(dx)) < toler
            break;
        end
    end
    elapsed = toc;
    Niter = k;
    
    delta = [0; x_state(1:nbus-1)];
    V = x_state(nbus:end);
end

% For the convenience of immediate reference, the helper functions are included underneath the body of the code

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
