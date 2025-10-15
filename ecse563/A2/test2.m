% IEEE 9 bus test system WECC representation
Sbase = 100; % MVA base

nfrom = [1 1 2]';
    
nto = [2 3 3]';

r = [0 0 0]';

x = [0.25 0.2 0.3]';

b = [0.25 0.2 0.3]';

% Power flow data

is = 1; ipq = [3]'; ipv = [1 2]'; V0 = [1 1]'; %[pv pv pq]

toler = 0.001; maxiter = 10;

Pd = [0.6 0 0.8]'; Qd = [0.2 0 0.3]'; % all known

Pg = [0 1 0]'; Qg = [0 0 0]'; %[? 1 0] [? ? 0] [pv pv pq]


%% part 1

% Number of buses - from the previous assignment
Y = admittance(nfrom, nto, r, x, b)

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