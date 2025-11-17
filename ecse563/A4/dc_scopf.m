
function out = dc_scopf(ifrom, ito, x, fmax, co, a, b, gmin, gmax, refbus, d, genbus, toler)
%DC_SCOPF DC OPF (intact-network security) with line limits via PTDFs.
%   Solves: min sum c0 + a*g + 0.5*b*g.^2
%   s.t.    sum(g) = sum(d), gmin<=g<=gmax, and |f|<=fmax (intact network)
%   using an active-set sequence that adds binding line constraints.
%   Returns dispatch, prices, flows, binding sets and costs.
%
%   out fields:
%     g, lambda_sys, mu_lines (map), flows, act.lines, LMP (all buses),
%     cost, cost_ED, cost_of_security, congestion_surplus
%
%   Reference: ECSE 563 notes (OPF + LMP, PTDF/SFT formulation).

n = max([ifrom(:); ito(:)]);
L = numel(ifrom);
G = numel(a);
a = a(:); b = b(:); co = co(:);
gmin = gmin(:); gmax = gmax(:); genbus = genbus(:);

% Build Bbus and PTDF (with refbus)
B = zeros(n,n);
bbr = 1./x(:);
for ell = 1:L
    i = ifrom(ell); j = ito(ell); b = bbr(ell);
    B(i,i) = B(i,i) + b; B(j,j) = B(j,j) + b;
    B(i,j) = B(i,j) - b; B(j,i) = B(j,i) - b;
end
mask = true(n,1); mask(refbus)=false;
X = zeros(n,n); X(mask,mask) = inv(B(mask,mask)); % maps injections to angles
% Incidence matrix (n x L): +1 at from, -1 at to
A = zeros(n,L); 
for ell=1:L, A(ifrom(ell),ell)=1; A(ito(ell),ell)=-1; end
C = diag(bbr) * A.';         % flow = C*delta
H = C*X;                     % PTDF: flow = H * injections
% Map generator outputs to bus injections P = M*g - d
M = zeros(n,G); for i=1:G, M(genbus(i), i) = 1; end
Hg = H * M;                  % (L x G)
fconst = H * (-d(:));        % (L x 1)

% Start from ED (no line limits)
Dtot = sum(d);
[tg, ~, lamED] = ed(co, a, b, gmin, gmax, Dtot, toler);
flows = Hg*tg + fconst;
% Active sets
act_lines = false(L,1);
% Work list: add most violated line until feasible
maxIter = 50; iter=0;
Q = diag(b); lin = a;   % quadratic/linear cost pieces
g = tg;
Aeq = [ones(1,G)]; beq = Dtot; names = {'balance'}; rows = [0];

while true
    iter = iter+1;
    if iter > maxIter, error('Active-set did not converge.'); end
    % Check line violations
    viol_up = flows - fmax(:);
    viol_lo = -fmax(:) - flows;
    [vmax_up, iup] = max(viol_up);
    [vmax_lo, ilo] = max(viol_lo);
    [vmax, which] = max([vmax_up, vmax_lo]);
    if vmax <= max(toler,1e-6)
        break; % feasible
    end
    if which==1
        ell = iup; sgn = +1; rhs = fmax(ell); tag = sprintf('line+%d',ell);
    else
        ell = ilo; sgn = -1; rhs = -fmax(ell); tag = sprintf('line-%d',ell);
    end
    % Add equality Hg(ell,:)*g = rhs - fconst(ell)
    Aeq = [Aeq; Hg(ell,:)];
    beq = [beq; rhs - fconst(ell)];
    names{end+1} = tag; rows(end+1)=ell;
    % Solve equality-constrained QP (unique g since #eq = #var may occur)
    KKT = [Q, -Aeq.'; Aeq, zeros(size(Aeq,1))];
    rhs_KKT = [-lin; beq];
    sol = KKT \ rhs_KKT;
    g = sol(1:G);
    y = sol(G+1:end);  % duals for [balance, added lines,...]
    flows = Hg*g + fconst;
    % Handle bound hits by pinning them and adding equalities
    hit_hi = find(g > gmax + 1e-8);
    hit_lo = find(g < gmin - 1e-8);
    for idx = hit_hi(:)'
        Aeq = [Aeq; unitrow(G, idx)];
        beq = [beq; gmax(idx)];
        names{end+1} = sprintf('gmax%d', idx); rows(end+1)=0;
    end
    for idx = hit_lo(:)'
        Aeq = [Aeq; unitrow(G, idx)];
        beq = [beq; gmin(idx)];
        names{end+1} = sprintf('gmin%d', idx); rows(end+1)=0;
    end
    % Re-solve after adding bounds if needed
    if ~isempty(hit_hi) || ~isempty(hit_lo)
        KKT = [Q, -Aeq.'; Aeq, zeros(size(Aeq,1))];
        rhs_KKT = [-lin; beq];
        sol = KKT \ rhs_KKT;
        g = sol(1:G);
        y = sol(G+1:end);
        flows = Hg*g + fconst;
    end
end

% Extract duals
lambda_sys = y(1);
% For line constraints, compute their multipliers in order of Aeq rows
mu_lines = zeros(L,1);
for k=2:numel(names)
    nm = names{k};
    if startsWith(nm,'line+')
        ell = rows(k); mu_lines(ell) = y(k);
    elseif startsWith(nm,'line-')
        ell = rows(k); mu_lines(ell) = y(k);
    end
end

% LMPs at all buses π = λ + sum_ell μ_ell * H(ell,:)
pi = lambda_sys + (mu_lines.' * H).';

% Congestion surplus: load payments minus gen payments
pay_load = sum(pi(:) .* d(:));
pay_gen  = sum(pi(genbus) .* g(:));
cong_surplus = pay_load - pay_gen;

C_ED = sum(co + a.*tg + 0.5*b.*(tg.^2));
C_OPF = sum(co + a.*g  + 0.5*b.*(g.^2));

out.g = g; out.lambda_sys = lambda_sys; out.mu_lines = mu_lines;
out.flows = flows; out.act = struct('lines', names);
out.LMP = pi; out.cost = C_OPF; out.cost_ED = C_ED;
out.cost_of_security = C_OPF - C_ED;
out.congestion_surplus = cong_surplus;
end

function e = unitrow(n, i)
e = zeros(1,n); e(i)=1;
end
