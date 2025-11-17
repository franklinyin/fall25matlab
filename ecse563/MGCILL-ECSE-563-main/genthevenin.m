function [Eeq, Zeq] = genthevenin(Y, Iint, id)

    % N is the number of rows / nodes
    N = size(Y,1);
    % reshape id into coulum vector - just in case :) 
    id = id(:);


    % Eq is Voc columns of id
    Voc = Y \ Iint;
    Eeq = Voc(id);

    % Zeq = Zii = ei' (Y^-1 ei)
    % id is a vector, so we can create a matrix of columns ei where i is
    % the index id

    k = size(id,1);
    Ei = sparse(N, k);
    for j = 1:k;
        Ei(id(j), j) = 1;
    end

    % Zeq = ei' Y^-1 ei 
    % we do it below in matrix Ei form Ei' Y^-1 Ei
    Zi = Y \ Ei;
    Zeq = Ei' * Zi;



end