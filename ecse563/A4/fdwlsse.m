% Q4 implementation
function [delta, V, Niter, elapsed] = fdwlsse(nfrom, nto, r, x, b, ...
    Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)
%FDWLSSE Weighted least-squares state estimation (AC) for small systems.
%   [delta, V, Niter, elapsed] = fdwlsse(nfrom, nto, r, x, b, ...
%       Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)
%   States:
%     delta : bus voltage angles (rad), size nb×1 (slack angle = 0)
%     V : bus voltage magnitudes (p.u.), size nb×1
%   Measurements:
%     Pinj : (mP×3) [bus, value(pu), sqrt(R)]
%     Qinj : (mQ×3) [bus, value(pu), sqrt(R)]
%     Pflow : (mPF×4)[from,to,value(pu),sqrt(R)]
%     Qflow : (mQF×4)[from,to,value(pu),sqrt(R)]
%     Vnode : (mV×3) [bus, value(pu), sqrt(R)]
%   Network:
%     nfrom, nto : line incidence (1-based bus indices)
%     r, x, b : line parameters (per unit); total line charging = j*b
%   Method:
%     Gauss-Newton WLS with numerical Jacobian and AC measurement model.
%     Slack bus is assumed to be bus 1 (angle fixed at 0).
%
%   Reference: ECSE 563 notes (Power Flow + State Estimation).

% Basic sizes
nbus = max(max(nfrom), max(nto));
nline = numel(nfrom);

% Build Ybus
Ybus = zeros(nbus);
yser = 1 ./ (r + 1j*x);
bsh = 1j * b / 2;

for l = 1:nline
    i = nfrom(l);
    k = nto(l);
    y = yser(l);
    bs = bsh(l);
    Ybus(i,i) = Ybus(i,i) + y + bs;
    Ybus(k,k) = Ybus(k,k) + y + bs;
    Ybus(i,k) = Ybus(i,k) - y;
    Ybus(k,i) = Ybus(k,i) - y;
end

% Extract measurements and build weight matrix
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

% Initial state guess
delta_state = zeros(nbus-1,1);
V = ones(nbus,1);
x = [delta_state; V];

% Iterative WLS
tic;
for k = 1:maxiter
    h = measurement_model(x, Ybus, nbus, nfrom, nto, yser, bsh, ...
                         Pinj, Qinj, Pflow, Qflow, Vnode);
    r = z - h;
    
    H = numerical_jacobian(x, h, Ybus, nbus, nfrom, nto, yser, bsh, ...
                          Pinj, Qinj, Pflow, Qflow, Vnode);
    
    G = H.' * (W * H);
    rhs = H.' * (W * r);
    
    dx = G \ rhs;
    
    x = x + dx;
    
    if max(abs(dx)) < toler
        break;
    end
end
elapsed = toc;
Niter = k;

% Unpack state
delta = [0; x(1:nbus-1)];
V = x(nbus:end);
end

% =====================================================================
function h = measurement_model(x, Ybus, nbus, nfrom, nto, yser, bsh, ...
                                Pinj, Qinj, Pflow, Qflow, Vnode)

delta_state = x(1:nbus-1);
V = x(nbus:end);

Va = [0; delta_state];
Vcmp = V .* exp(1j*Va);

Iinj = Ybus * Vcmp;
Sbus = Vcmp .* conj(Iinj);
Pbus = real(Sbus);
Qbus = imag(Sbus);

h = [];

% P injections
if ~isempty(Pinj)
    for k = 1:size(Pinj,1)
        i = Pinj(k,1);
        h = [h; Pbus(i)];
    end
end

% Q injections
if ~isempty(Qinj)
    for k = 1:size(Qinj,1)
        i = Qinj(k,1);
        h = [h; Qbus(i)];
    end
end

% P flows
if ~isempty(Pflow)
    for k = 1:size(Pflow,1)
        i = Pflow(k,1);
        j = Pflow(k,2);
        idx = find(nfrom==i & nto==j, 1);
        if isempty(idx)
            idx = find(nfrom==j & nto==i, 1);
            i2 = j; j2 = i;
        else
            i2 = i; j2 = j;
        end
        y = yser(idx);
        bs = bsh(idx);
        Iij = (Vcmp(i2) - Vcmp(j2)) * y + Vcmp(i2) * bs;
        Sij = Vcmp(i2) * conj(Iij);
        h = [h; real(Sij)];
    end
end

% Q flows
if ~isempty(Qflow)
    for k = 1:size(Qflow,1)
        i = Qflow(k,1);
        j = Qflow(k,2);
        idx = find(nfrom==i & nto==j, 1);
        if isempty(idx)
            idx = find(nfrom==j & nto==i, 1);
            i2 = j; j2 = i;
        else
            i2 = i; j2 = j;
        end
        y = yser(idx);
        bs = bsh(idx);
        Iij = (Vcmp(i2) - Vcmp(j2)) * y + Vcmp(i2) * bs;
        Sij = Vcmp(i2) * conj(Iij);
        h = [h; imag(Sij)];
    end
end

% V magnitudes
if ~isempty(Vnode)
    for k = 1:size(Vnode,1)
        i = Vnode(k,1);
        h = [h; abs(Vcmp(i))];
    end
end

end

% =====================================================================
function H = numerical_jacobian(x, h0, Ybus, nbus, nfrom, nto, yser, bsh, ...
                                 Pinj, Qinj, Pflow, Qflow, Vnode)

nstate = numel(x);
m = numel(h0);
H = zeros(m, nstate);
eps_fd = 1e-6;

for j = 1:nstate
    xpert = x;
    xpert(j) = xpert(j) + eps_fd;
    hpert = measurement_model(xpert, Ybus, nbus, nfrom, nto, yser, bsh, ...
                              Pinj, Qinj, Pflow, Qflow, Vnode);
    H(:,j) = (hpert - h0) / eps_fd;
end
end
