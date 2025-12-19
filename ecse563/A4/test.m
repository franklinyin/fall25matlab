% a4 run file
clear; clc;

% data for q1&q2
c0 = [800; 400; 400];
a  = [20; 23; 56];
b  = [0.030; 0.035; 0.040];
gmin = [10; 10; 10];
gmax = [250; 300; 270];
toler = 5e-1;  % MW

loads = [300, 400, 600, 700]; % q2

%% q1
fprintf('Q1: ED\n');
for d = loads
    [g_ed,C_ed,lam_ed] = ed(c0,a,b,gmin,gmax,d,toler);
    fprintf('  d=%4.0f MW -> lambda=%.3f $/MWh, g=[%.2f %.2f %.2f] MW, C=%.2f $/h\n', ...
        d, lam_ed, g_ed, C_ed);
end

% production visualization
q1_visualization(loads, c0, a, b, gmin, gmax, toler);

% fixed cost recovery analysis
q1_min_output(loads, c0, a, b, gmin, gmax, toler);


%% q2
fprintf('\nQ2: UC\n');
for d = loads
    [u_uc,g_uc,C_uc,lam_uc] = uc(c0,a,b,gmin,gmax,d,toler);
    prof_uc = lam_uc*g_uc - (c0.*u_uc + a.*g_uc + 0.5*b.*g_uc.^2);
    fprintf('- d=%4.0f: u=%s, lambda=%.3f, g=[%.2f,%.2f,%.2f],\n    C=%.2f, total profit=%.2f\n\n', ...
        d, mat2str(u_uc.'), lam_uc, g_uc, C_uc, sum(prof_uc));
end

%% q3
A4Q3_scopf_data;      % expects variables as provided
out = dc_scopf(ifrom, ito, x, fmax, d, co, a, b, gmin, gmax, ngen, is);
fprintf('\nQ3: DC SCOPF with intact line limits\n');
fprintf('  g* = [%g %g %g] MW\n', out.g);
fprintf('  cost = %.2f $/h\n', out.C);
fprintf('  Merchandizing surplus = %.2f $/h\n', out.MS);
disp('  LMPs by bus:');
disp(out.LMP.');

%% ---------- Problem 4 data (A4Q4_wlsse_data.m) ----------
A4Q4_wlsse_data;
[delta, V, Niter, tsec] = fdwlsse(nfrom, nto, r, x, b, Pinj, Qinj, Pflow, Qflow, Vnode, 1e-4, 50);
fprintf('\nQ4: WLS State Estimation\n');
fprintf('  Iterations=%d, time=%.4fs\n', Niter, tsec);
fprintf('  Angles (rad):\n'); disp(delta.');
fprintf('  Voltages (pu):\n'); disp(V.');