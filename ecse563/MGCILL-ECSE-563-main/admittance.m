function Y = admittance(nfrom, nto, r, x, b)

    % N is the number of rows / nodes
    N = max(max(nfrom), max(nto));
    % Number of branches
    M = length(nfrom);

    % Make the A matrix
    % size of matrix is N x M (nodes x branches/lines)
    Id = eye(N);
    A = Id(1:N, nfrom) - Id(1:N, nto);
    Asize = size(A);

    % What to do with A if 2 columns are identical --> Parralel lines

    % Calculate admittance per branch. 
    % Make sure to use ./ for element wise divison
    zl = r + 1j*x;
    yl = 1 ./ zl;
    Yb = diag(yl);

    % Include shunt
    Y = A*Yb*A'+ 1j*0.5*diag(abs(A)*b);



end




