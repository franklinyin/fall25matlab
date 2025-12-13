function out = dc_scopf(ifrom, ito, x, fmax, d, co, a, b, gmin, gmax, ngen, is)
% DC SCOPF for the modified IEEE 9-bus system (intact network)
%
% Inputs:
%   ifrom, ito  : line "from" and "to" bus indices (nlines x 1)
%   x           : line reactances (p.u. on Sbase) (nlines x 1)
%   fmax        : line MW flow limits (nlines x 1)
%   d           : bus demands (MW) (nbus x 1)
%   co,a,b      : generator cost coefficients
%                 C_i(g_i) = co_i + a_i g_i + 0.5 b_i g_i^2
%   gmin,gmax   : generator min/max (MW) (ng x 1)
%   ngen        : bus index of each generator (ng x 1)
%   is          : reference (slack) bus index
%
% Outputs in struct out:
%   g       : generator outputs (MW)
%   delta   : bus angles (rad), ref bus = 0
%   f       : line flows (MW), positive ifrom -> ito
%   C       : total generation cost ($/h)
%   LMP     : locational marginal prices ($/MWh) per bus
%   MS      : merchandizing (congestion) surplus ($/h)
%   lambda  : lambda struct from quadprog (KKT multipliers)
%   pinj    : net injections per bus (MW)
%   ifrom, ito, refbus : echoed inputs

%% Basic dimensions
nbus   = max([ifrom; ito]);
nlines = length(ifrom);
ng     = length(co);
genbus = ngen(:);   % bus index per generator

%% Build line susceptances and Bbus
bline = 1 ./ x;     % susceptance per line (1/x)

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
keep   = setdiff(1:nbus, refbus);   % all buses except reference
Bred   = B(keep, keep);

%% Incidence matrix A (line -> bus)
A = zeros(nlines, nbus);
for ell = 1:nlines
    A(ell, ifrom(ell)) =  1;
    A(ell, ito(ell))   = -1;
end


% Line flow matrix: f = F * delta_red
F = diag(bline) * A(:, keep);   % nlines x (nbus-1)

%% Generator incidence matrix G (bus -> gen)
G = zeros(nbus, ng);
for k = 1:ng
    G(genbus(k), k) = 1;
end

Gk  = G(keep, :);     % rows for non-ref buses
dk  = d(keep);

%% QP matrices: minimize 0.5 z' H z + fvec' z
% z = [g; delta_red]

H    = blkdiag(diag(b), zeros(length(keep)));    % angles have no cost
fvec = [a(:); zeros(length(keep),1)];

% Nodal balances (for non-ref buses):
%   Bred * delta_red = Gk * g - dk
% → [-Gk  Bred] [g; delta_red] = -dk
Aeq_nodal = [-Gk, Bred];
beq_nodal = -dk;

% Global power balance: sum_i g_i = sum_n d_n
Aeq_bal = [ones(1,ng), zeros(1,length(keep))];
beq_bal = sum(d);

% Full equality set
Aeq = [Aeq_nodal;
       Aeq_bal];
beq = [beq_nodal;
       beq_bal];

% Line limits: -fmax <= F * delta_red <= fmax
Aineq = [zeros(nlines, ng),  F;
         zeros(nlines, ng), -F];
bineq = [fmax(:); fmax(:)];

% Variable bounds
lb = [gmin(:);  -inf(length(keep),1)];
ub = [gmax(:);   inf(length(keep),1)];

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

delta = zeros(nbus,1);
delta(keep) = delta_red;           % reference bus angle = 0

% Line flows
f = F * delta_red;                 % MW, oriented ifrom -> ito

% Cost
C = sum(co + a(:).*g + 0.5*b(:).*g.^2);

%% LMPs (dual variables of nodal balance equations)
lambda_bus = zeros(nbus,1);

% In lambda.eqlin, rows correspond to rows of Aeq:
% first (nbus-1) are nodal balances, last one is global balance
lambda_nodal = lambda.eqlin(1:length(keep));

lambda_bus(keep) = lambda_nodal;

% set ref bus price as average of others (or any consistent value)
lambda_bus(refbus) = mean(lambda_bus(keep));

LMP = lambda_bus;

%% Congestion / merchandizing surplus
% Net injections: generation positive, load positive
pinj = G * g - d;   % nbus x 1

% MS = - sum_i lambda_i * p_i
MS = -sum(LMP .* pinj);

%% Pack outputs
out.g       = g;
out.delta   = delta;
out.f       = f;
out.C       = C;
out.LMP     = LMP;
out.MS      = MS;
out.lambda  = lambda;
out.pinj    = pinj;
out.ifrom   = ifrom;
out.ito     = ito;
out.refbus  = refbus;

end
