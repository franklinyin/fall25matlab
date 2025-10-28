% Load IEEE 9-bus system data
run('ieee9_A2.m');


%% Q1

% Number of buses - from the previous assignment
Y = admittance(nfrom, nto, r, x, b);

[V_nrpf, delta_nrpf, Psl_nrpf, Qgv_nrpf, N_nrpf, time_nrpf] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Q1: \n');
fprintf('Converged in %d iterations and took %f seconds.\n', N_nrpf, time_nrpf);
disp('Voltage magnitudes (p.u.):');
disp(V_nrpf);
disp('Voltage angles (radians):');
disp(delta_nrpf);
disp('Active power at slack bus (MW):');
disp(Psl_nrpf);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv_nrpf);


%% Q2

fprintf('Q2: \n');
[V_decpf, delta_decpf, Psl_decpf, Qgv_decpf, N_decpf, time_decpf] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N_decpf, time_decpf);
disp('Voltage magnitudes (p.u.):');
disp(V_decpf);
disp('Voltage angles (radians):');
disp(delta_decpf);
disp('Active power at slack bus (MW):');
disp(Psl_decpf);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv_decpf);


%% Q3

fprintf('Q3: \n');
[V_fastdecpf, delta_fastdecpf, Psl_fastdecpf, Qgv_fastdecpf, N_fastdecpf, time_fastdecpf] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N_fastdecpf, time_fastdecpf);
disp('Voltage magnitudes (p.u.):');
disp(V_fastdecpf);
disp('Voltage angles (radians):');
disp(delta_fastdecpf);
disp('Active power at slack bus (MW):');
disp(Psl_fastdecpf);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv_fastdecpf);


%% Q4


fprintf('Q4: \n');
[delta_dcpf, Psl_dcpf, Pf_dcpf] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase);

% Display the results
disp('Voltage angles (radians):');
disp(delta_dcpf);
disp('Active power at slack bus (MW):');
disp(Psl_dcpf);
disp('Active power at transmission line (MW):');
disp(Pf_dcpf);



%% Q5 Compare AC vs DC active line flows

Sij_complex = compute_acpf(V_fastdecpf, delta_fastdecpf, Y, nfrom, nto, Sbase);

fprintf('Q5: Compare AC vs DC active line flows\n');
fprintf('\nLine  from-to    P_AC(MW)     S_AC(MVA)      P_DC(MW)     |P_AC|(MW)    |S_AC|(MVA)\n');
for e = 1:numel(nfrom)
    Pf_ac = real(Sij_complex(e));
    Sf_ac = abs(Sij_complex(e));
    Pf_diff = Pf_ac - Pf_dcpf(e);
    Sf_diff = Sf_ac - abs(Pf_dcpf(e));
    
    fprintf('%2d    %d  -> %d   %9.3f     %9.3f     %9.3f     %9.3f     %9.3f\n', e, nfrom(e), nto(e), Pf_ac, Sf_ac, Pf_dcpf(e), abs(Pf_ac), abs(Sf_ac));
end


%% Q6 Feasibility (at bus 7)

run('ieee9_A2.m');

fprintf('Q6: Feasibility (at bus 7)\n');

P7_vals = 0:20:400; Q7_vals = -200:20:200;
feas = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, P7_vals, Q7_vals, toler, maxiter);
fprintf('\nFeasibility map computed on %dx%d grid. Feasible points: %d\n', numel(P7_vals), numel(Q7_vals), nnz(feas));

figure; % plot
imagesc(P7_vals, Q7_vals, double(feas));
set(gca,'YDir','normal'); axis tight
colormap([0.9 0.5 0.5; 0.5 0.9 0.5]);
cb = colorbar('Ticks',[0 1],'TickLabels',{'infeasible','feasible'});
xlabel('P_{d7} (MW)'); ylabel('Q_{d7} (MVAr)');
title('Feasibility region for bus 7 (varying load at node 7)'); grid on

% Save the plot as PNG
saveas(gcf, 'feasibility_map_bus7.png');


%% Q7 Security analysis (N-1)
fprintf('Q7: Security analysis (N-1)\n');
cont = [4 5; 4 9; 5 6; 6 7; 7 8; 8 9];
out = pf_security_analysis(nfrom, nto, r, x, b, Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, cont, toler, maxiter);
fprintf('\nN-1 contingencies (acceptable range 0.95–1.05 p.u.):\n');
for k = 1:numel(out)
    fprintf('Outage %d-%d: Vmin=%.3f, Vmax=%.3f\n', out(k).pair(1), out(k).pair(2), out(k).Vmin, out(k).Vmax);
end