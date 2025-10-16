% Q6
function feas = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, P7_vals, Q7_vals, toler, maxiter)
%PF_FEASIBILITY_MAP  Boolean grid (|P7| x |Q7|) true if NRPF converges.
%   sweeps Pd(7) and Qd(7) and tests NR convergence. (No Q-limits enforced.)
n = numel(Pd); feas = false(numel(P7_vals), numel(Q7_vals));
for a = 1:numel(P7_vals)
    for b = 1:numel(Q7_vals)
        Pd_try = Pd; Qd_try = Qd;
        Pd_try(7) = P7_vals(a);
        Qd_try(7) = Q7_vals(b);
        [~, ~, ~, ~, N, ~] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd_try, Qd_try, V0, Sbase, toler, maxiter);
        feas(a,b) = (N < maxiter);   % crude convergence test
    end
end
end
