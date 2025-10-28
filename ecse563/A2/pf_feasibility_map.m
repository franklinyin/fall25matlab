% Q6 helper function
function feas = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, P7_vals, Q7_vals, toler, maxiter)
    % pf_feasibility_map - compute feasibility map for varying load at bus 7
    % Input:
    %   Y - admittance matrix
    %   is - slack bus index
    %   ipq, ipv - PQ and PV bus indices
    %   Pg, Qg - generation vectors
    %   Pd, Qd - demand vectors
    %   V0 - initial voltage magnitudes
    %   Sbase - base power (MVA)
    %   P7_vals, Q7_vals - load variation ranges for bus 7
    %   toler, maxiter - convergence parameters
    % Output:
    %   feas - feasibility matrix (boolean)

    feas = false(numel(P7_vals), numel(Q7_vals));
    for a = 1:numel(P7_vals)
        for b = 1:numel(Q7_vals)
            Pd_try = Pd; Qd_try = Qd;
            Pd_try(7) = P7_vals(a);
            Qd_try(7) = Q7_vals(b);
            [V, ~, ~, ~, N, ~] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd_try, Qd_try, V0, Sbase, toler, maxiter);
            if N < maxiter
                if all(abs(V)>=0.95 & abs(V)<=1.05)
                    feas(a,b) = true;
                end
            end
        end
    end
end
