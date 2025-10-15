%% Q5
function [IT, VNF] = genfault(YN, YF, IintN, IintF, idN, idF)
    % genfault - Perform generalized fault calculation
    % Input:
    %   YN - admittance matrix of the healthy network
    %   YF - admittance matrix of the fault network
    %   IintN - pre-fault internal currents in the healthy network
    %   IintF - pre-fault internal currents in the fault network
    %   idN - node indices in the healthy network
    %   idF - node indices in the fault network
    % Output:
    %   IT - tie line currents
    %   VNF - node voltages in the healthy network with the fault network connected
    
    k = size(idN,1);

    permF = [idF', setdiff(1:size(YF, 1), idF)]'; % permutation
    permN = [idN', setdiff(1:size(YN,1), idN)]';

    YFp = YF(permF,permF);    IFp = IintF(permF); % dense, < 25 lines but readable

    YF11 = YFp(1:k,1:k);     YF12 = YFp(1:k,k+1:end); 
    YF21 = YFp(k+1:end,1:k); YF22 = YFp(k+1:end,k+1:end);

    IF1 = IFp(1:k);  IF2 = IFp(k+1:end);

    X = YF22 \ [YF21, IF2]; % YF22^-1*YF21 and YF22^-1*IF2 zipped, then unpack below:
    YFeq = YF11 - YF12*X(:,1:k); % (10) p.41 in Fault Analysis lhs
    IFeq = IF1 - YF12*X(:,k+1); % (10) p.41 in Fault Analysis rhs

    YNp = YN(permN,permN); INp = IintN(permN);

    YN11 = YNp(1:k,1:k);     YN12 = YNp(1:k,k+1:end);
    YN21 = YNp(k+1:end,1:k); YN22 = YNp(k+1:end,k+1:end);

    rhs = [INp(1:k)+IFeq; INp(k+1:end)];
    VHp = [YN11+YFeq, YN12; YN21, YN22] \ rhs; % Rearranging (11) p. 43 in Fault Analysis

    % recover
    VNF = zeros(size(YN,1),1); 
    VNF(permN) = VHp;

    V1 = VHp(1:k);
    IT = YFeq*V1 - IFeq; % last equation on p. 49 in Fault Analysis
end