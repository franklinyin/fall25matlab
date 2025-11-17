function [delta, Psl, Pf, time] = dcpf(nfrom, nto, x, is, Pg, Pd, Sbase)
tic;

% Initialization
nbus = max([nfrom; nto]);
nline = length(nfrom);

Pg = Pg(:);
Pd = Pd(:);
x  = x(:);
nfrom = nfrom(:);
nto   = nto(:);

%  Build susceptance matrix B' (imaginary part only, no shunts) ---
B = sparse(nbus, nbus);

for k = 1:nline
    i = nfrom(k); 
    j = nto(k);
    bij = 1 / x(k);
    B(i,i) = B(i,i) + bij;
    B(j,j) = B(j,j) + bij;
    B(i,j) = B(i,j) - bij;
    B(j,i) = B(j,i) - bij;
end

%  Convert active power to per-unit
Pinj = (Pg - Pd) / Sbase;

%  Solve linear system - remove slack bus row/column
idx = setdiff((1:nbus)', is);
delta = zeros(nbus,1);
delta(idx) = B(idx, idx) \ Pinj(idx);

% Compute slack power injection from solved angles
Pcalc = B * delta;
Psl   = Pcalc(is) * Sbase;     % in MW

% Compute line active power flows (MW), from nfrom --> nto
Pf = ((delta(nfrom) - delta(nto)) ./ x) * Sbase;

time = toc;

end
