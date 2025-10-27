% Load IEEE 9-bus system data
run('ieee9_A2.m');


%% part 1

% Number of buses - from the previous assignment
Y = admittance(nfrom, nto, r, x, b);

% Call the Newton-Raphson power flow function
[V_nrpf, delta_nrpf, Psl_nrpf, Qgv_nrpf, N_nrpf, time_nrpf] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N_nrpf, time_nrpf);
disp('Voltage magnitudes (p.u.):');
disp(V_nrpf);
disp('Voltage angles (radians):');
disp(delta_nrpf);
disp('Active power at slack bus (MW):');
disp(Psl_nrpf);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv_nrpf);



%% part 2


% Call the Newton-Raphson power flow function
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


%% part 3


% Call the Newton-Raphson power flow function
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


%% part 4


% Call the Newton-Raphson power flow function
[delta_dcpf, Psl_dcpf, Pf_dcpf] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase);

% Display the results
disp('Voltage angles (radians):');
disp(delta_dcpf);
disp('Active power at slack bus (MW):');
disp(Psl_dcpf);
disp('Active power at transmission line (MW):');
disp(Pf_dcpf);


% Q7
% === Security analysis (N-1) ===
cont = [4 5; 4 9; 5 6; 6 7; 7 8; 8 9];
out = pf_security_analysis(nfrom, nto, r, x, b, Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, cont, toler, maxiter);
fprintf('\nN-1 contingencies (acceptable range 0.95–1.05 p.u.):\n');
for k = 1:numel(out)
    fprintf('Outage %d-%d: Vmin=%.3f, Vmax=%.3f\n', out(k).pair(1), out(k).pair(2), out(k).Vmin, out(k).Vmax);
end