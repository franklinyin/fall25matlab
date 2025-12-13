
% RUN_A4  — quick validation harness for Assignment 4
clear; clc;

% ---------- Problem 1 & 2 data (Table 1) ----------
c0 = [800; 400; 400];
a  = [20; 23; 56];
b  = [0.030; 0.035; 0.040];
gmin = [10; 10; 10];
gmax = [250; 300; 270];
toler = 5e-1;  % MW

loads = [300, 400, 600, 700];
fprintf('Problem 1: ED results\n');
for d = loads
    [g,C,lam] = ed(c0,a,b,gmin,gmax,d,toler);
    fprintf('  d=%4.0f MW -> lambda=%.3f $/MWh, g=[%.2f %.2f %.2f] MW, C=%.2f $/h\n', ...
        d, lam, g, C);
end

fprintf('\nProblem 2: UC (static, full enumeration)\n');
for d = loads
    [u,g,C,lam] = uc(c0,a,b,gmin,gmax,d,toler);
    prof = lam*g - (c0.*u + a.*g + 0.5*b.*g.^2);
    fprintf('  d=%4.0f: u=%s, lambda=%.3f, g=[%.2f %.2f %.2f], C=%.2f, total profit=%.2f\n', ...
        d, mat2str(u.'), lam, g, C, sum(prof));
end

% ---------- Problem 3 data (A4Q3_scopf_data.m) ----------
A4Q3_scopf_data;      % expects variables as provided
out = dc_scopf(ifrom, ito, x, fmax, d, co, a, b, gmin, gmax, ngen, is);
fprintf('\nProblem 3: DC SCOPF with intact line limits\n');
fprintf('  g* = [%g %g %g] MW\n', out.g);
fprintf('  cost = %.2f $/h\n', out.C);
fprintf('  Merchandizing surplus = %.2f $/h\n', out.MS);
disp('  LMPs by bus:');
disp(out.LMP.');

% ---------- Problem 4 data (A4Q4_wlsse_data.m) ----------
A4Q4_wlsse_data;
[delta, V, Niter, tsec] = fdwlsse(nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, 1e-4, 50);
fprintf('\nProblem 4: WLS State Estimation\n');
fprintf('  Iterations=%d, time=%.4fs\n', Niter, tsec);
fprintf('  Angles (rad):\n'); disp(delta.');
fprintf('  Voltages (pu):\n'); disp(V.');
