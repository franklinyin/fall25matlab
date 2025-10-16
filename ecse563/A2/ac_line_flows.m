% Q5
function [Pf, Qf, Pt, Qt, S_from, S_to] = ac_line_flows(nfrom, nto, r, x, b, V, delta, Sbase)
%AC_LINE_FLOWS  Branch powers using series admittance + line charging b/2 per end.
%   Returns P,Q from 'from' to 'to' (Pf,Qf) and reverse (Pt,Qt), plus |S| magnitudes.
m = numel(nfrom);
Pf = zeros(m,1); Qf = zeros(m,1);
Pt = zeros(m,1); Qt = zeros(m,1);
S_from = zeros(m,1); S_to = zeros(m,1);
for e = 1:m
    i = nfrom(e); k = nto(e);
    z   = r(e) + 1i*x(e);
    y   = 1/z;
    bsh = 1i*b(e)/2;
    Vi  = V(i)*exp(1i*delta(i));
    Vk  = V(k)*exp(1i*delta(k));
    Iik = (Vi - Vk)*y + Vi*bsh;   % current leaving i towards k
    Iki = (Vk - Vi)*y + Vk*bsh;   % current leaving k towards i
    Sij = Vi*conj(Iik)*Sbase;     % MVA
    Sji = Vk*conj(Iki)*Sbase;
    Pf(e) = real(Sij); Qf(e) = imag(Sij);
    Pt(e) = real(Sji); Qt(e) = imag(Sji);
    S_from(e) = abs(Sij); S_to(e) = abs(Sji);
end
end
