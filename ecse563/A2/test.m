% IEEE 9 bus test system WECC representation
Sbase = 100; % MVA base

nfrom = [1 4 5 3 6 7 8 8 9]';
    
nto = [4 5 6 6 7 8 2 9 4]';

r = [0 0.017 0.039 0 0.0119 0.0085 0 0.032 0.01]';

x = [0.0576 0.092 0.17 0.0586 0.1008 0.072 0.0625 0.161 0.085]';

b = [0 0.158 0.358 0 0.209 0.149 0 0.306 0.176]';

% Power flow data

is = 1; ipq = [4 5 6 7 8 9]'; ipv = [1 2 3]'; V0 = [1 1 1]';

toler = 0.001; maxiter = 10;

Pd = [0 0 0 0 90 0 100 0 125]'; Qd = [0 0 0 0 30 0 35 0 50]';

Pg = [0 163 85 0 0 0 0 0 0]'; Qg = [0 0 0 0 0 0 0 0 0]';


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



%% Q5

% Set the parameters for the Newton-Raphson power flow
Y = admittance(nfrom, nto, r, x, b); % Build admittance matrix from system data
V0 = ones(9, 1); % Starting voltages at 1 p.u. for all buses
toler = 0.001; % Convergence tolerance
maxiter = 20; % Maximum iterations

% Define the inputs for Pg, Pd, Qg, Qd as previously defined
[~, ~, Psl_NR, ~, ~, ~] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);

% Calculate line flows in AC power flow
Pac = zeros(length(nfrom), 1);
for i = 1:length(nfrom)
    Pac(i) = real(Sbase * Y(nfrom(i), nto(i)) * (V0(nfrom(i)) - V0(nto(i))) * conj(V0(nfrom(i))));
end

%step 2
fprintf('Line\tFrom\tTo\tPdc (MW)\tPac (MW)\n');
for i = 1:length(nfrom)
    fprintf('%d\t%d\t%d\t%.2f\t%.2f\n', i, nfrom(i), nto(i), Pf(i), Pac(i));
end


%% part 6

% Assuming previous setup of the IEEE 9 bus system
% Define ranges for Pd7 and Qd7
Pd7_range = linspace(0, 300, 100);  % Pd7 from 0 to 300 MW
Qd7_range = linspace(-100, 100, 100);  % Qd7 from -100 to 100 MVar

% Preallocate for results
feasible = zeros(length(Pd7_range), length(Qd7_range));

% Load original system data
Pg_original = [163; 85; 0; 0; 0; 0; 0; 0; 0];  % Original generation in MW
Pd_original = [0; 0; 0; 0; 90; 0; 100; 0; 125]; % Original demand in MW
Qd_original = [0; 0; 0; 0; 30; 0; 35; 0; 50];  % Original reactive demand in MVar

for i = 1:length(Pd7_range)
    for j = 1:length(Qd7_range)
        Pd = Pd_original;
        Qd = Qd_original;
        Pd(7) = Pd7_range(i);
        Qd(7) = Qd7_range(j);
        
        try
            % Run Newton-Raphson power flow
            [~, ~, Psl, ~, ~, ~] = nrpf(Y, is, ipq, ipv, Pg_original, Qg, Pd, Qd, V0, Sbase, toler, maxiter);
            feasible(i, j) = 1;  % Mark this point as feasible
        catch
            feasible(i, j) = 0;  % Mark this point as infeasible if NR fails
        end
    end
end

% Plotting the feasibility region
figure;
imagesc(Pd7_range, Qd7_range, feasible');
xlabel('Active Power Demand at Bus 7 (MW)');
ylabel('Reactive Power Demand at Bus 7 (MVar)');
title('Feasibility Region for Load at Bus 7');
colorbar;
axis xy;


%% Part 7
% List of contingencies
contingencies = [4 5; 4 9; 5 6; 6 7; 7 8; 8 9];

% Initialize system parameters
% Assuming Y, Pg, Pd, Qg, Qd, V0, Sbase have been set up correctly

% Store original admittance matrix
Y_original = Y;

% Store results
results = struct();

for i = 1:size(contingencies, 1)
    % Modify the admittance matrix for the contingency
    Y = Y_original;
    Y(contingencies(i,1), contingencies(i,2)) = 0;
    Y(contingencies(i,2), contingencies(i,1)) = 0;
    Y(contingencies(i,1), contingencies(i,1)) = sum(Y_original(contingencies(i,1), :));
    Y(contingencies(i,2), contingencies(i,2)) = sum(Y_original(contingencies(i,2), :));
    
    % Solve power flow
    [V, ~, ~, ~, ~, ~] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, 0.001, 20);
    
    % Record the results
    results(i).line = contingencies(i,:);
    results(i).voltages = V;
    results(i).acceptable = all(V >= 0.95 & V <= 1.05);
end

% Display results
for i = 1:size(contingencies, 1)
    fprintf('Contingency Line %d-%d:\n', results(i).line(1), results(i).line(2));
    disp(results(i).voltages');
    if results(i).acceptable
        fprintf('Voltages are within the acceptable range.\n\n');
    else
        fprintf('Voltages are out of the acceptable range!\n\n');
    end
end
