%% Q2
function Z = impedance(nfrom, nto, r, x, b)
    % impedance - calculate the impedance matrix Z
    % Input:
    %   nfrom - vector of sending buses
    %   nto - vector of receiving buses
    %   r, x, b - vectors of resistance, reactance and shunt susceptance
    % Output:
    %   Z - impedance matrix

    Y = admittance(nfrom, nto, r, x, b); % get admittance matrix from Q1
    Z = linsolve(Y, eye(size(Y))); % inverse admittance matrix using identity matrix
end