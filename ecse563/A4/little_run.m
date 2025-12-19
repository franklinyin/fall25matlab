clear
% ECSE 563 assigmenet 4
% Ali Seifeldin
% https://github.com/Bakalala/MGCILL-ECSE-563

% Question 1

clear; clc; close all;

% Load table 1 data.
run("A4Q3_scopf_data.m")

%% Load levels and tolerance
loads = [300 400 600 700];    % MW
toler = 0.5;                  

nGen   = numel(co);
nCases = numel(loads);

G          = zeros(nGen, nCases);   % generator outputs (MW)
C_total    = zeros(1, nCases);      % total cost ($/h)
Lambda     = zeros(1, nCases);      % system lambda ($/MWh)
Profit     = zeros(nGen, nCases);   % per-unit profit ($/h)
TotalProfit = zeros(1, nCases);     % total profit ($/h)




%% Run ED for each load level and store results
for k = 1:nCases
    d = loads(k);

    [gk, Ck, lambdak] = ed(co, a, b, gmin, gmax, d, toler);

    G(:, k)       = gk;
    C_total(k)    = Ck;
    Lambda(k)     = lambdak;

    Cost_i       = co + a .* gk + 0.5 * b .* (gk.^2);  % each unit's cost
    Revenue_i    = lambdak * gk;                       % λ * g_i
    Profit_i     = Revenue_i - Cost_i;                 % profit_i

    Profit(:,k)   = Profit_i;
    TotalProfit(k) = sum(Profit_i);
end

%% Build a single results table (with profits)
Results = table( ...
    loads(:), ...        % Load
    Lambda(:), ...       % Lambda
    G(1,:).', ...        % g1
    G(2,:).', ...        % g2
    G(3,:).', ...        % g3
    Profit(1,:).', ...   % Profit1
    Profit(2,:).', ...   % Profit2
    Profit(3,:).', ...   % Profit3
    TotalProfit(:), ...  % Total system profit
    C_total(:), ...      % Total cost
    'VariableNames', { ...
        'Load_MW','Lambda', ...
        'g1_MW','g2_MW','g3_MW', ...
        'Profit1','Profit2','Profit3', ...
        'TotalProfit','TotalCost'} ...
);

%% Display the table
fprintf('=== ED algorithm Table ===\n');
disp(Results)


%% Plot incremental cost curves and operating points for one load
% Choose which load to illustrate graphically:
d_plot  = 600;                            % MW (change if you like)
k_plot  = find(loads == d_plot, 1)

for k = 1:nCases

    d_plot  = loads(k);                            % MW (change if you like)

    if ~isempty(k)
        lambda_plot = Lambda(k);
        g_plot      = G(:, k);
    
        figure; hold on; grid on;
        title(sprintf('Incremental cost curves at d = %d MW (lambda = %.3f)', ...
            d_plot, lambda_plot));
        xlabel('Generator output g_i (MW)');
        ylabel('Incremental cost dC_i/dg_i ($/MWh)');
    
        for i = 1:nGen
            g_range = linspace(gmin(i), gmax(i), 200);
            mc      = a(i) + b(i) .* g_range;   % dC_i/dg_i = a_i + b_i*g_i
    
            plot(g_range, mc, 'DisplayName', sprintf('Unit %d MC', i));
            plot(g_plot(i), lambda_plot, 'o', 'HandleVisibility', 'off');
            text(g_plot(i), lambda_plot, sprintf('  U%d', i), ...
                'VerticalAlignment', 'bottom');
        end
    
        yline(lambda_plot, '--', 'DisplayName', '\lambda');
        legend('Location', 'best');
    end
end

%% Compute minimum outputs to recover fixed costs at each lambda
% Solve for each unit and each lambda:
%   lambda * g = c0 + a*g + 0.5*b*g^2
% -> 0.5*b*g^2 + (a - lambda)*g + c0 = 0

g_break_even = NaN(nGen, nCases);

for k = 1:nCases
    lam = Lambda(k);

    for i = 1:nGen

        % Solve A*g^2 + B*g + C = 0 using MATLAB's roots()
        coeffs = [0.5*b(i),  a(i) - lam,  co(i)];
        r = roots(coeffs);

        % Select real, positive roots only
        r = r(imag(r)==0 & r > 0);

        if ~isempty(r)
            g_break_even(i, k) = min(r);   % smallest positive real root
        end
    end
end


%% Table: Minimum outputs to recover fixed costs at each lambda
% (solution of lambda*g = c0 + a*g + 0.5*b*g^2)

BreakEven = table( ...
    loads(:), ...             % Load for each case
    Lambda(:), ...            % Lambda for each case
    g_break_even(1,:).', ...  % g_be1 for each case (column)
    g_break_even(2,:).', ...  % g_be2 for each case
    g_break_even(3,:).', ...  % g_be3 for each case
    'VariableNames', {'Load_MW','Lambda','g_be1_MW','g_be2_MW','g_be3_MW'} ...
);

%'Minimum outputs to recover fixed costs at each lambda'
% solution of lambda*g = c0 + a*g + 0.5*b*g^2

fprintf('=== ED algorithm BREAKEVEN Table ===\n');
disp(BreakEven)

% The table shows that for 300MW, we never break even as the prices are too
% low
% For 700 MW, we breakeven with all 3 generators.

fprintf('=== ED algorithm Table ===\n');
disp(Results)

% Question 2

clear; clc; close all;

% Load table 1 data.
run("A4Q3_scopf_data.m")

loads = [300 400 600 700];    % MW

toler = 0.5;  

nGen   = numel(co);
nCases = numel(loads);

U       = zeros(nGen, nCases);   % generator scheduling status
G       = zeros(nGen, nCases);   % generator outputs (MW)
C_total = zeros(1, nCases);      % total cost ($/h)
Lambda  = zeros(1, nCases);      % system lambda ($/MWh)
Profit  = zeros(nGen, nCases);   % Profit (MW)

for k = 1:nCases
    d = loads(k);
    [u,g,C,lambda] = uc(co,a,b,gmin,gmax,d,toler);

    % Per-unit costs & profits at that lambda:
    Cost_i   = u .* co + a .* g + 0.5 * b .* (g.^2);
    Revenue  = lambda * g;           % λ * g_i
    Profit_i = Revenue - Cost_i;

    U(:,k)      = u;
    G(:,k)      = g;
    C_total(k)  = C;
    Lambda(k)   = lambda;
    Profit(:,k) = Profit_i;

end

fprintf('=== UC algorithm Table ===\n');
TotalProfit = sum(Profit, 1);  % Calculate total profit for each case

%% Build a single results table
Results = table( ...
    loads(:), ...         % Load
    Lambda(:), ...        % Lambda
    U(1,:).', ...         % U1
    U(2,:).', ...         % U2
    U(3,:).', ...         % U3
    G(1,:).', ...         % g1
    G(2,:).', ...         % g2
    G(3,:).', ...         % g3
    Profit(1,:).', ...    % Profit1
    Profit(2,:).', ...    % Profit2
    Profit(3,:).', ...    % Profit3
    TotalProfit(:), ...   % <-- ADD THIS
    C_total(:), ...       % Total cost
    'VariableNames', { ...
        'Load_MW','Lambda', ...
        'U1','U2','U3', ...
        'g1_MW','g2_MW','g3_MW', ...
        'Profit1','Profit2','Profit3', ...
        'TotalProfit', ...
        'TotalCost'} ...
);

%% Display the table
fprintf('=== UC algorithm Table ===\n');
disp(Results)


% Using the UC algorithm, we can see that at the 400MW load, we make a
% profit, as opposed to losing money using the ED algo. This is due to the
% fact that we completely shut off generator 3, as opposed to keeping it
% on. 
% We get a sense of that using from the breakeven calculations we did, where it
% is impossible to have any of the generator breakeven until we reach at least
% 300MW. Instead of keeping it on we shut it off.
% we still lose money when we generate 300 MW but we minimized the losses
% from ~ -1100 to ~ 400 per hour
% 600 and 700 MW loads were profitable in both solutions. since all
% generators need to be so the solutions are identical
% we can notice that Generator 3 for both solutions with 600 MW is
% unprofitable, even if the system as a whole is profitable. 

% Question 3

clear; clc; close all;

% Load table 1 data.
run("A4Q3_scopf_data.m")

%% A4Q3 main script


%% 1) Economic dispatch ignoring network (from Q1)
Dt = sum(d);  % total demand (MW)
[g_ed, C_ed, lambda_sys] = ed(co,a,b,gmin,gmax,Dt,dp);

%% 2) DC "SCOPF" with network
out = dc_scopf(ifrom, ito, x, fmax, ...
                     d, co, a, b, gmin, gmax, ngen, is);

g_scopf = out.g;
C_scopf = out.C;
LMP     = out.LMP;
MS      = out.MS;
flows   = out.f;

%% 3) Cost of security
cost_of_security = C_scopf - C_ed;

%% 4) Display key results
disp('--- Economic Dispatch (no network) ---');
disp(table((1:length(g_ed))', g_ed, ...
    'VariableNames', {'Gen','P_ED_MW'}));
fprintf('Total ED cost = %.2f $/h\n', C_ed);

disp(' ');
disp('--- DC SCOPF (with 9-bus network) ---');
disp(table((1:length(g_scopf))', g_scopf, ...
    'VariableNames', {'Gen','P_SCOPF_MW'}));
fprintf('Total SCOPF cost = %.2f $/h\n', C_scopf);
fprintf('Cost of security = %.2f $/h\n\n', cost_of_security);

disp('Bus LMPs ($/MWh):');
disp(table((1:length(LMP))', LMP, ...
    'VariableNames', {'Bus','LMP'}));

disp('Line flows (MW):');
disp(table((1:length(flows))', ifrom, ito, flows, fmax, ...
    'VariableNames', {'Line','From','To','Flow','Fmax'}));

fprintf('Merchandizing surplus (congestion surplus) = %.2f $/h\n', MS);

% The SCOPF solution redispatches generators to ensure that all line-flow 
% constraints are satisfied under the DC network model. Compared to the ED case, 
% generator 2 reaches its maximum output because it is located electrically 
% closer to load centers, reducing power transfer across heavily loaded lines. 
% Generator 1 decreases slightly, and generator 3 increases output to relieve congestion 
% on lines such as line 9–4, which reaches its limit of 105 MW. 
% These adjustments increase total cost relative to ED because the system 
% must operate in a less economical configuration to obey transmission limits.

% The cost of security of approximately 585 $/h represents the economic penalty 
% of enforcing transmission constraints. In the ED case, power flows freely, 
% but in SCOPF several lines become congested—particularly line 9–4 and lines connected to buses 7 and 8. 
% Congestion prevents the system from using the cheapest generation pattern and 
% forces redispatch toward higher-cost units. 

% The LMP values vary significantly across the network, ranging from near 0 $/MWh at bus 4 
% to over 68 $/MWh at bus 9. These variations reflect congestion patterns: buses electrically 
% downstream of congested lines face higher marginal costs because supplying additional 
% load there would further violate line limits. Buses near cheap generation with unconstrained 
% transfer capability (such as bus 4) have low LMPs. The merchandizing surplus of 12,421 $/h arises 
% from LMP differences across the network and represents congestion rents collected by the system operator. 
% This surplus is a direct consequence of the binding line constraints observed in the SCOPF solution.





% Question 4

clear; clc; close all;

% Load table 1 data.
run("A4Q4_wlsse_data.m")

%% Run state estimator
[delta, V, Niter, elapsed] = fdwlsse( ...
    nfrom, nto, r, x, b, ...
    Pinj, Qinj, Pflow, Qflow, Vnode, ...
    toler, maxiter);

%% Display results
fprintf('=== WLS State Estimation Results ===\n');
fprintf('Iterations: %d   Time: %.6f s\n\n', Niter, elapsed);

nbus = numel(V);
fprintf('Bus   Voltage(pu)   Angle(deg)\n');
for k = 1:nbus
    fprintf('%3d   %10.4f   %10.4f\n', k, V(k), delta(k)*180/pi);
end

% Comparingg the output to table 4. 
% V1, V3 and V4 are respectivelty 1,1,1.05. 
% Our state estimated 1.0062,  1.0101, 1.0582
% This is reaspnabe, with less than a 1% difference
% we can also infer that bus 3 has a generator with a + angle
% we also notice bus 4 has a negative angle bus has a generator
% that is because the load is higher than what is generated so that is ok


% Calculation for observability

% 9 sates to observe, 4 angles and 5 voltages.
% Angle at slack is 0

N = max(max(nfrom), max(nto));
Id = eye(N);
A = Id(1:N, nfrom) - Id(1:N, nto);
A = A'

% For Active power

Maa = [1 1 0 0 0 0 0;
    0 0 -1 1 1 0 0;
    0 -1 0 -1 0 1 -1;
    1 0 0 0 0 0 0]

%Remove slack node
A_noslack = A(:, 2:end)
Haa = Maa * A_noslack

Gaa = transpose(Haa) * Haa
rankGaa = rank(Gaa)

% Rank is 4 --> all 4 angles are observable

% For reactive power

Mrr = [1 1 0 0 0 0 0;
    0 0 -1 1 1 0 0;
    0 -1 0 -1 0 1 -1;
    1 0 0 0 0 0 0 ;]

% add measurement for V2
Hrr = [Mrr * A;
        0 1 0 0 0]

Hrr = transpose(Hrr) * Hrr
rankGrr = rank(Hrr)

totalRank = rankGaa + rankGrr
% Rank is 5 --> all 5 magnitudes are observable

%Rank 4 + 5 = 9 --> fully obserable. 
% If we remove 1 measurement we would
%reduce in rank for sure, so we would not have full observability






