function [delta, Psl, Pf] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase)
    % dcpf - DC power flow (angles only; line MW flows from nfrom->nto)
    % Input:
    %   nfrom, nto - line connection vectors
    %   x - line reactance vector
    %   is - slack bus index
    %   Pg, Pd - generation and demand vectors
    %   Sbase - base power (MVA)
    % Output:
    %   delta - voltage angles (radians)
    %   Psl - slack bus active power (MW)
    %   Pf - line active power flows (MW)
%
% Assumptions: |V|≈1 p.u., small angles, R≈0, shunts ignored.
% Angles solve  B' * delta = Pinj   (slack removed), and line flows are
%   Pf_l = (delta_i - delta_k)/x_l * Sbase   in MW. :contentReference[oaicite:19]{index=19} :contentReference[oaicite:20]{index=20}
%
n = max(max(nfrom), max(nto));
Bbus = zeros(n);
m = numel(nfrom);
for e = 1:m
    i = nfrom(e); k = nto(e);
    bij = 1/x(e);
    Bbus(i,i) = Bbus(i,i) + bij;
    Bbus(k,k) = Bbus(k,k) + bij;
    Bbus(i,k) = Bbus(i,k) - bij;
    Bbus(k,i) = Bbus(k,i) - bij;
end

Pinj = (Pg - Pd)/Sbase;     % p.u.
keep = setdiff(1:n, is);
Bred = Bbus(keep, keep);
Pred = Pinj(keep);

dred = Bred \ Pred;
delta = zeros(n,1);
delta(keep) = dred;   % slack angle = 0

% Slack active power (MW) — by balance (no losses in DC):
Psl = Pd(is) - sum( Pg(keep) - Pd(keep) );

% Line flows (MW) from nfrom->nto
Pf = zeros(m,1);
for e = 1:m
    i = nfrom(e); k = nto(e);
    Pf(e) = ( delta(i) - delta(k) ) / x(e) * Sbase;
end
end
