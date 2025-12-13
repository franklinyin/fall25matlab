function [delta, V, Niter, elapsed] = fdwlsse( ...
    nfrom, nto, r, x, b, ...
    Pinj, Qinj, Pflow, Qflow, Vnode, ...
    toler, maxiter)
%FDWLSSE  Weighted least-squares state estimation (AC) for small systems.
%
%   [delta, V, Niter, elapsed] = fdwlsse(nfrom, nto, r, x, b, ...
%       Pinj, Qinj, Pflow, Qflow, Vnode, toler, maxiter)
%
%   States:
%       delta : bus voltage angles (rad), size nb×1 (slack angle = 0)
%       V     : bus voltage magnitudes (p.u.), size nb×1
%
%   Measurements:
%       Pinj  : (mP×3) [bus, value(pu), sqrt(R)]
%       Qinj  : (mQ×3) [bus, value(pu), sqrt(R)]
%       Pflow : (mPF×4)[from,to,value(pu),sqrt(R)]
%       Qflow : (mQF×4)[from,to,value(pu),sqrt(R)]
%       Vnode : (mV×3) [bus, value(pu), sqrt(R)]
%
%   Network:
%       nfrom, nto : line incidence (1-based bus indices)
%       r, x, b    : line parameters (per unit); total line charging = j*b
%
%   Method:
%       Gauss-Newton WLS with numerical Jacobian and AC measurement model.
%       Slack bus is assumed to be bus 1 (angle fixed at 0).

% ---------- basic sizes ----------
nbus  = max(max(nfrom), max(nto));
nline = numel(nfrom);

% ---------- build Ybus ----------
Ybus = zeros(nbus);
yser = 1 ./ (r + 1j*x);
bsh  = 1j * b / 2;

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

% J(x) = (z - h(x))^T W (z - h(x))
% stack to form this objective function

% ---------- stack measurements ----------
z      = [];
bus1   = [];      % first bus index per measurement
bus2   = [];      % second bus (for flows)
mtype  = [];      % 1=Pinj, 2=Qinj, 3=Pflow, 4=Qflow, 5=V
sigma2 = [];      % variances

% P injections
if ~isempty(Pinj)
    z      = [z ; Pinj(:,2)];
    bus1   = [bus1 ; Pinj(:,1)];
    bus2   = [bus2 ; zeros(size(Pinj,1),1)];
    mtype  = [mtype; 1*ones(size(Pinj,1),1)];
    sigma2 = [sigma2 ; Pinj(:,3).^2];
end

% Q injections
if ~isempty(Qinj)
    z      = [z ; Qinj(:,2)];
    bus1   = [bus1 ; Qinj(:,1)];
    bus2   = [bus2 ; zeros(size(Qinj,1),1)];
    mtype  = [mtype; 2*ones(size(Qinj,1),1)];
    sigma2 = [sigma2 ; Qinj(:,3).^2];
end

% P flows
if ~isempty(Pflow)
    z      = [z ; Pflow(:,3)];
    bus1   = [bus1 ; Pflow(:,1)];
    bus2   = [bus2 ; Pflow(:,2)];
    mtype  = [mtype; 3*ones(size(Pflow,1),1)];
    sigma2 = [sigma2 ; Pflow(:,4).^2];
end

% Q flows
if ~isempty(Qflow)
    z      = [z ; Qflow(:,3)];
    bus1   = [bus1 ; Qflow(:,1)];
    bus2   = [bus2 ; Qflow(:,2)];
    mtype  = [mtype; 4*ones(size(Qflow,1),1)];
    sigma2 = [sigma2 ; Qflow(:,4).^2];
end

% V magnitudes
if ~isempty(Vnode)
    z      = [z ; Vnode(:,2)];
    bus1   = [bus1 ; Vnode(:,1)];
    bus2   = [bus2 ; zeros(size(Vnode,1),1)];
    mtype  = [mtype; 5*ones(size(Vnode,1),1)];
    sigma2 = [sigma2 ; Vnode(:,3).^2];
end

m = numel(z);                    % number of measurements
W = diag(1 ./ sigma2);           % weighting matrix

% ---------- initial state guess ----------
delta_state = zeros(nbus-1,1);   % angles except slack
V           = ones(nbus,1);      % magnitudes
x           = [delta_state; V];  % stacked state vector
nstate      = numel(x);

% function handle for measurement model
measfun = @(xx) measurement_model(xx, Ybus, nbus, ...
                                  nfrom, nto, yser, bsh, ...
                                  mtype, bus1, bus2);

% ---------- iterative WLS ----------
tic;
for k = 1:maxiter
    h = measfun(x);
    r = z - h;                   % residual

    % numerical Jacobian
    H = numerical_jacobian(measfun, x, h);

    G   = H.' * (W * H);         % gain matrix
    rhs = H.' * (W * r);

    dx = G \ rhs;

    x  = x + dx;

    if max(abs(dx)) < toler
        break;
    end
end
elapsed = toc;
Niter   = k;

% ---------- unpack state ----------
delta = [0; x(1:nbus-1)];
V     = x(nbus:end);

end

% =====================================================================
function h = measurement_model(x, Ybus, nbus, ...
                               nfrom, nto, yser, bsh, ...
                               mtype, bus1, bus2)
% build predicted measurements h(x)
% yser vector of series admittance
% bsh vecotre of shunt admittance

delta_state = x(1:nbus-1);
V           = x(nbus:end);

Va = [0; delta_state];                  % include slack angle
Vcmp = V .* exp(1j*Va);                 % complex voltages

% bus injections
Iinj = Ybus * Vcmp;
Sbus = Vcmp .* conj(Iinj);              % complex power injections
Pbus = real(Sbus);
Qbus = imag(Sbus);

h = zeros(numel(mtype),1);

for k = 1:numel(mtype)
    t = mtype(k);
    i = bus1(k);
    j = bus2(k);

    switch t
        case 1      % Pinj at bus i
            h(k) = Pbus(i);
        case 2      % Qinj at bus i
            h(k) = Qbus(i);
        case 3      % Pflow from i to j
            idx = find(nfrom==i & nto==j, 1);
            if isempty(idx)
                idx = find(nfrom==j & nto==i, 1);
                i2  = j; j  = i;
            else
                i2 = i;
            end
            y  = yser(idx);
            bs = bsh(idx);
            Iij = (Vcmp(i2) - Vcmp(j)) * y + Vcmp(i2) * bs;
            Sij = Vcmp(i2) * conj(Iij);
            h(k) = real(Sij);
        case 4      % Qflow from i to j
            idx = find(nfrom==i & nto==j, 1);
            if isempty(idx)
                idx = find(nfrom==j & nto==i, 1);
                i2  = j; j  = i;
            else
                i2 = i;
            end
            y  = yser(idx);
            bs = bsh(idx);
            Iij = (Vcmp(i2) - Vcmp(j)) * y + Vcmp(i2) * bs;
            Sij = Vcmp(i2) * conj(Iij);
            h(k) = imag(Sij);
        case 5      % V magnitude at bus i
            h(k) = abs(Vcmp(i));
    end
end
end

% =====================================================================
function H = numerical_jacobian(measfun, x, h0)
% simple finite-difference Jacobian
nstate = numel(x);
m      = numel(h0);
H      = zeros(m, nstate);
eps_fd = 1e-6;

for j = 1:nstate
    xpert    = x;
    xpert(j) = xpert(j) + eps_fd;
    hpert    = measfun(xpert);
    H(:,j)   = (hpert - h0) / eps_fd;
end
end
