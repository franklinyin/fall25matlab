% Load IEEE 9-bus system data
run('ieee9_A2.m');


%% part 1

% Number of buses - from the previous assignment
Y = admittance(nfrom, nto, r, x, b);

% Call the Newton-Raphson power flow function
[V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N, time);
disp('Voltage magnitudes (p.u.):');
disp(V);
disp('Voltage angles (radians):');
disp(delta);
disp('Active power at slack bus (MW):');
disp(Psl);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv);



%% part 2


% Call the Newton-Raphson power flow function
[V, delta, Psl, Qgv, N, time] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N, time);
disp('Voltage magnitudes (p.u.):');
disp(V);
disp('Voltage angles (radians):');
disp(delta);
disp('Active power at slack bus (MW):');
disp(Psl);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv);


%% part 3


% Call the Newton-Raphson power flow function
[V, delta, Psl, Qgv, N, time] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Display the results
fprintf('Converged in %d iterations and took %f seconds.\n', N, time);
disp('Voltage magnitudes (p.u.):');
disp(V);
disp('Voltage angles (radians):');
disp(delta);
disp('Active power at slack bus (MW):');
disp(Psl);
disp('Reactive power at PV buses (Mvar):');
disp(Qgv);


%% part 4


% Call the Newton-Raphson power flow function
[delta, Psl, Pf] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase);

% Display the results
disp('Voltage angles (radians):');
disp(delta);
disp('Active power at slack bus (MW):');
disp(Psl);
disp('Active power at transmission line (MW):');
disp(Pf);
