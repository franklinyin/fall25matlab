% Q3 implementation
function out = dc_scopf(ifrom, ito, x, fmax, d, co, a, b, gmin, gmax, ngen, is)
%DC_SCOPF DC security-constrained optimal power flow (intact network)
%   Solves: min sum c0_i + a_i*g_i + 0.5*b_i*g_i^2
%   s.t.    nodal power balance, generation limits, and line flow limits
%   Inputs:
%     ifrom, ito : line "from" and "to" bus indices (nlines x 1)
%     x : line reactances (p.u.) (nlines x 1)
%     fmax : line MW flow limits (nlines x 1)
%     d : bus demands (MW) (nbus x 1)
%     co,a,b : generator cost coefficients
%     gmin,gmax : generator min/max (MW) (ng x 1)
%     ngen : bus index of each generator (ng x 1)
%     is : reference (slack) bus index
%   Outputs (struct out):
%     g : generator outputs (MW)
%     delta : bus angles (rad), ref bus = 0
%     f : line flows (MW), positive ifrom -> ito
%     C : total generation cost ($/h)
%     LMP : locational marginal prices ($/MWh) per bus
%     MS : merchandizing (congestion) surplus ($/h)
%     lambda : lambda struct from quadprog (KKT multipliers)
%     pinj : net injections per bus (MW)
%     ifrom, ito, refbus : echoed inputs
%
%   Reference: ECSE 563 notes (OPF + LMP, DC approximation).

%% Basic dimensions
nbus   = max([ifrom; ito]);
nlines = length(ifrom);
ng     = length(co);
genbus = ngen(:);

%% Build line susceptances and Bbus
bline = 1 ./ x;

B = zeros(nbus);
for ell = 1:nlines
    i = ifrom(ell);
    j = ito(ell);
    B(i,i) = B(i,i) + bline(ell);
    B(j,j) = B(j,j) + bline(ell);
    B(i,j) = B(i,j) - bline(ell);
    B(j,i) = B(j,i) - bline(ell);
end

% Reduced B (remove reference bus)
refbus = is;
keep   = setdiff(1:nbus, refbus);
Bred   = B(keep, keep);

%% Incidence matrix A (line -> bus)
A = zeros(nlines, nbus);
for ell = 1:nlines
    A(ell, ifrom(ell)) = 1;
    A(ell, ito(ell))   = -1;
end

% Line flow matrix: f = F * delta_red
F = diag(bline) * A(:, keep);

%% Generator incidence matrix G (bus -> gen)
G = zeros(nbus, ng);
for k = 1:ng
    G(genbus(k), k) = 1;
end

Gk = G(keep, :);
dk = d(keep);

%% QP matrices: minimize 0.5 z' H z + fvec' z
% z = [g; delta_red]
Hgen = diag(b(:));
H    = blkdiag(Hgen, zeros(length(keep)));
fvec = [a(:); zeros(length(keep),1)];

%% Equality constraints

% Nodal balances (for non-ref buses):
%   Bred * delta_red = Gk * g - dk
% -> [-Gk  Bred] [g; delta_red] = -dk
Aeq_nodal = [-Gk, Bred];
beq_nodal = -dk;

% Global power balance: sum_i g_i = sum_n d_n
Aeq_bal = [ones(1,ng), zeros(1,length(keep))];
beq_bal = sum(d);

% Full equality set
Aeq = [Aeq_nodal; Aeq_bal];
beq = [beq_nodal; beq_bal];

%% Line limits: -fmax <= F * delta_red <= fmax
Aineq = [zeros(nlines, ng),  F;
         zeros(nlines, ng), -F];
bineq = [fmax(:); fmax(:)];

%% Variable bounds
lb = [gmin(:); -inf(length(keep),1)];
ub = [gmax(:);  inf(length(keep),1)];

%% Solve QP with quadprog
options = optimoptions('quadprog','Display','off');

[z, ~, exitflag, ~, lambda] = quadprog(H, fvec, Aineq, bineq, ...
                                       Aeq, beq, lb, ub, [], options);

if exitflag <= 0
    error('quadprog did not converge in dc_scopf (exitflag = %d)', exitflag);
end

%% Extract solution
g         = z(1:ng);
delta_red = z(ng+1:end);

delta        = zeros(nbus,1);
delta(keep)  = delta_red;
% refbus angle is zero by construction

% Line flows
f = F * delta_red;

% Cost
C = sum(co + a(:).*g + 0.5*b(:).*g.^2);

%% LMPs (shadow prices of bus power balance equations)

% lambda.eqlin corresponds to rows of Aeq:
%  1..(nbus-1): nodal balances for buses in 'keep'
%  last       : global power balance
lambda_eqlin  = lambda.eqlin;
nkeep         = length(keep);
lambda_nodalM = lambda_eqlin(1:nkeep);      % for nodal balances (non-ref buses)
lambda_balM   = lambda_eqlin(nkeep + 1);    % for global balance

% Sensitivity of optimal cost w.r.t. load at each bus (slide 11):
%  For non-ref bus i in 'keep': d_i appears with -1 in its nodal balance
%  and +1 in the global balance RHS.
%  Using quadprog sign convention, this gives:
%     LMP_i = lambda_nodalM(i) - lambda_balM
%  For the reference bus, load appears only in the global balance:
%     LMP_ref = -lambda_balM
LMP = zeros(nbus,1);
LMP(keep)  = lambda_nodalM - lambda_balM;
LMP(refbus) = -lambda_balM;

%% Congestion / merchandizing surplus

% Net injections: generation positive, load positive
pinj = G * g - d;

% MS = - sum_i LMP_i * p_i
MS = -sum(LMP .* pinj);

%% Pack outputs
out.g      = g;
out.delta  = delta;
out.f      = f;
out.C      = C;
out.LMP    = LMP;
out.MS     = MS;
out.lambda = lambda;
out.pinj   = pinj;
out.ifrom  = ifrom;
out.ito    = ito;
out.refbus = refbus;
end
