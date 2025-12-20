% Observability Proof

N = max(max(nfrom), max(nto));
Id = eye(N);
A = (Id(1:N, nfrom) - Id(1:N, nto))';

% active power observability
Maa = [1 1 0 0 0 0 0; 0 0 -1 1 1 0 0; 0 -1 0 -1 0 1 -1; 1 0 0 0 0 0 0];
Haa = Maa * A(:, 2:end);  % Remove slack node
Gaa = Haa' * Haa;
rankGaa = rank(Gaa);  % Rank 4: all 4 angles observable

% reactive power observability
Mrr = [1 1 0 0 0 0 0; 0 0 -1 1 1 0 0; 0 -1 0 -1 0 1 -1; 1 0 0 0 0 0 0];
Hrr = [Mrr * A; 0 1 0 0 0];  % Add V2 measurement
Grr = Hrr' * Hrr;
rankGrr = rank(Grr);  % Rank 5: all 5 magnitudes observable

totalRank = rankGaa + rankGrr;  % Rank 9: fully observable
