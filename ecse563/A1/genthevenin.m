%% Q4
function [Eeq, Zeq] = genthevenin(Y, Iint, id)
    % genthevenin - calculate generalized Thevenin equivalent
    % Input:
    %   Y - admittance matrix of the network
    %   Iint - vector of pre-fault internal currents
    %   id - vector of node indices for Thevenin equivalent
    % Output:
    %   Eeq - Thevenin equivalent voltages
    %   Zeq - Thevenin equivalent impedance matrix

    % reorder the Y matrix and Iint
    nNode = length(id);
    newOrderY = [id, setdiff(1:size(Y, 1), id)];
    Y_permuted = Y(newOrderY, newOrderY);
    Iint_permuted = Iint(newOrderY);

    % subdivid the matrix into 4 parts, the Iint into 2 parts
    YF11 = Y_permuted(1:nNode, 1:nNode);
    YF12 = Y_permuted(1:nNode, nNode+1:end);
    YF21 = Y_permuted(nNode+1:end, 1:nNode);
    YF22 = Y_permuted(nNode+1:end, nNode+1:end);
    IF1 = Iint_permuted(1:nNode);
    IF2 = Iint_permuted(nNode+1:end);

    % solving for YFeq and IFeq
    YFeq = YF11 - YF12 * linsolve(YF22,YF21);
    IFeq = IF1 - YF12 * linsolve(YF22, IF2);

    % solving for Zeq and Eeq
    Zeq = linsolve(YFeq,eye(nNode));
    Eeq = linsolve(YFeq,IFeq);
end