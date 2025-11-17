function Z = impedance(nfrom, nto, r, x, b)


    Y = admittance(nfrom, nto, r, x, b);   % N x N

    % Solve for Z since YZ = Identity matrix
    N = size(Y,1);
    Z = Y \ eye(N);



end




