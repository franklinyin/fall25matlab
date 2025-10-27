%% Q1 from assignment 1
function Y = admittance(nfrom, nto, r, x, b)
    % admittance - calculate the admittance matrix Y
    % Input:
    %   nfrom - vector of sending buses
    %   nto - vector of receiving buses
    %   r, x, b - vectors of resistance, reactance and shunt susceptance
    % Output:
    %   Y - admittance matrix

    n = max(max(nfrom), max(nto)); % determine matrix size using the biggest node number
    Y = zeros(n, n); % initialize admittance matrix
    z = r + 1i*x; % impedance of each line: z = r + j*x
    y = 1 ./ z; % admittance of each line y = 1/z
    
    % Fill admittance matrix
    for k = 1:size(nfrom,1)
        Y(nfrom(k), nto(k)) = Y(nfrom(k), nto(k)) - y(k); % Off-diagonal
        Y(nto(k), nfrom(k)) = Y(nto(k), nfrom(k)) - y(k); % Symmetry
        Y(nfrom(k), nfrom(k)) = Y(nfrom(k), nfrom(k)) + y(k) + 1i*b(k)/2; % self admittance
        Y(nto(k), nto(k)) = Y(nto(k), nto(k)) + y(k) + 1i*b(k)/2;
    end
end