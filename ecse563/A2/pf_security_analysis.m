% Q7
function out = pf_security_analysis(nfrom, nto, r, x, b, baseY, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, contingencies, toler, maxiter)
%PF_SECURITY_ANALYSIS  For each (i,k) outage in 'contingencies', remove line and solve NR.
%   Returns struct with fields: pair, converged, N, Vmin, Vmax.
out = struct('pair',{},'converged',{},'N',{},'Vmin',{},'Vmax',{});
for c = 1:size(contingencies,1)
    i = contingencies(c,1); k = contingencies(c,2);
    % remove line i-k (either orientation)
    mask = true(numel(nfrom),1);
    for e = 1:numel(nfrom)
        if (nfrom(e)==i && nto(e)==k) || (nfrom(e)==k && nto(e)==i), mask(e)=false; end
    end
    Yc = admittance(nfrom(mask), nto(mask), r(mask), x(mask), b(mask));
    [V, delta, ~, ~, N, ~] = nrpf(Yc, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);
    out(c).pair      = [i k];
    out(c).converged = (N < maxiter);
    out(c).N         = N;
    out(c).Vmin      = min(V);
    out(c).Vmax      = max(V);
end
end
