function [If, Vf] = fault(Y, Iint, idfault, Zf)

    % N is the number of rows / nodes
    N = size(Y,1);

    % Create ei vector, with 1 at idfault
    ei = sparse(N,1); 
    ei(idfault) = 1;

    % Calculate Voc and Zi
    Voc = Y \ Iint;
    Zi = Y \ ei;

    % Calculate Vf and If as per formulas
    Vf = Voc - Zi * (ei' * Voc) / (ei' * Zi + Zf);
    If = (ei' * Voc) / (ei' * Zi + Zf);



end




