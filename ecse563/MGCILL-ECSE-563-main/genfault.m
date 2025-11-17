function [IT, VNF] = genfault(YN, YF, IintN, IintF, idN, idF)
    
    % Re-index Y and Iint based on id.
    [YN, IintN] = arrange_from_id(YN, IintN, idN);
    [YF, IintF] = arrange_from_id(YF, IintF, idF);


    % number of nodes connected between healthy and fault network. 
    % Need this to 'break' the matrix into 4 subsection
    m = size(idN,1);

    % Breakdown the Y matrix and I into Y11/Y12/Y21/Y22 and I1/I2 for
    % healthy network
    YN11 = YN(1:m,1:m);     
    YN12 = YN(1:m,m+1:end);
    YN21 = YN(m+1:end,1:m); 
    YN22 = YN(m+1:end,m+1:end);
    IN1  = IintN(1:m);      
    IN2  = IintN(m+1:end);

    % Breakdown the Y matrix and I into Y11/Y12/Y21/Y22 and I1/I2 for fault
    % network
    YF11 = YF(1:m,1:m);     
    YF12 = YF(1:m,m+1:end);
    YF21 = YF(m+1:end,1:m); 
    YF22 = YF(m+1:end,m+1:end);
    IF1  = IintF(1:m);      
    IF2  = IintF(m+1:end);


    % Find Y_eq and I_eq for fault network
    YF_eq = YF11 - YF12 * (YF22 \ YF21);
    IF_eq = IF1  - YF12 * (YF22 \ IF2);

    % YN equal to = [YN11+ YF_eq, YN12 ; YN21, YN22]
    % IN is equal to = [IN1 + IF_eq; IN2]
    VNF = [YN11+ YF_eq, YN12 ; YN21, YN22] \ [IN1 + IF_eq; IN2];

    % Take the first m rows to get V1 to calculate IT into the 
    V1 = VNF(1:m);
    % IT = Yeq V1 - Ieq - current flowing into healthy network
    IT = YF_eq*V1 - IF_eq;


              


end


% Function to reaarange matrix to 'bring up' terms in id
function [Y,Iint] = arrange_from_id(Y, Iint, id)
    
    % number of rows
    N = size(Y,1);
    
    % Get all indexes 1 ... N, remove values in id and add to the front of
    % the vector.
    all_id = (1:N).';
    all_id = setdiff(all_id,id);
    all_id = [id;all_id];

    % Re-index Y and Iint based on the new indexing order.
    Y = Y(all_id,all_id);
    Iint = Iint(all_id);


end