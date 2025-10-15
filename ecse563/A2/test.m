clear; clc;

% === Load network ===
ieee9_A2;                        % your provided data file (variables in workspace)
Y = admittance(nfrom, nto, r, x, b);

% Force assignment settings
toler   = 1e-4;
maxiter = 20;

% PV magnitudes are given for ipv order in the file:
VO = V0;

fprintf('\n=== Newton–Raphson (NRPF) ===\n');
[V, delta, Psl, Qgv, NitNR, tNR] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, VO, Sbase, toler, maxiter);
fprintf('Converged in %d iters (%.4g s). Slack P = %.3f MW\n', NitNR, tNR, Psl);
pv_wo_slack = setdiff(ipv, is);
for k = 1:numel(pv_wo_slack)
    i = pv_wo_slack(k); fprintf('  Qg@bus %d = %.3f MVAr\n', i, Qgv(i));
end
fprintf('Voltage magnitudes (p.u.):\n'); disp([(1:numel(V))' V]);
fprintf('Voltage angles (deg):\n');       disp([(1:numel(delta))' rad2deg(delta)]);

% Quick software validation: mismatches below toler at PQ/PV, and power balance
[Pf, Qf, Pt, Qt, Sfrom, Sto] = ac_line_flows(nfrom, nto, r, x, b, V, delta, Sbase);
Ploss = sum(Pf + Pt);
assert(Ploss > 0);  % positive losses
Pg_tot = sum(Pg) + Psl;
Pd_tot = sum(Pd);
assert(abs(Pg_tot - Pd_tot - Ploss) < 1e-3, 'Power balance check failed.');

fprintf('\n=== Decoupled PF (DPF) ===\n');
[Vd, deltad, Psld, Qgvd, NitD, tD] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, VO, Sbase, toler, maxiter);
fprintf('Converged in %d iters (%.4g s). Slack P = %.3f MW\n', NitD, tD, Psld);

fprintf('\n=== Fast-Decoupled PF (FDPF) ===\n');
[Vf, deltaf, Pslf, Qgvf, NitF, tF] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, VO, Sbase, toler, maxiter);
fprintf('Converged in %d iters (%.4g s). Slack P = %.3f MW\n', NitF, tF, Pslf);

fprintf('\n=== DC PF (DCPF) ===\n');
[delta_dc, Psl_dc, Pf_dc] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase);
fprintf('Slack P (DC) = %.3f MW\n', Psl_dc);

% Compare AC vs DC active line flows
fprintf('\nLine  from-to    P_AC(MW)     |S|_from(MVA)   P_DC(MW)\n');
for e = 1:numel(nfrom)
    fprintf('%2d    %d  -> %d   %9.3f     %9.3f     %9.3f\n', e, nfrom(e), nto(e), Pf(e), Sfrom(e), Pf_dc(e));
end

% === Feasibility (bus 7) ===
P7_vals = 0:20:400; Q7_vals = -400:20:400;
feas = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, VO, Sbase, P7_vals, Q7_vals, toler, maxiter);
fprintf('\nFeasibility map computed on %dx%d grid. Feasible points: %d\n', numel(P7_vals), numel(Q7_vals), nnz(feas));

% === Security analysis (N-1) ===
cont = [4 5; 4 9; 5 6; 6 7; 7 8; 8 9];
out = pf_security_analysis(nfrom, nto, r, x, b, Y, is, ipq, ipv, Pg, Qg, Pd, Qd, VO, Sbase, cont, toler, maxiter);
fprintf('\nN-1 contingencies (acceptable range 0.95–1.05 p.u.):\n');
for k = 1:numel(out)
    fprintf('Outage %d-%d: Vmin=%.3f, Vmax=%.3f\n', out(k).pair(1), out(k).pair(2), out(k).Vmin, out(k).Vmax);
end
