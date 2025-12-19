% Observability Proof

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






