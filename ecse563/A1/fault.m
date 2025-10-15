%% Q3
function [If, Vf] = fault(Y, Iint, idfault, Zf)
    % fault - calculate fault current and node voltages
    % Input:
    %   Y - admittance matrix of the network
    %   Iint - vector of pre-fault internal currents
    %   idfault - index of the faulted node
    %   Zf - fault impedance
    % Output:
    %   If - fault current
    %   Vf - node voltages during fault

    % setting the e_i
    N = length(Iint); 
    e_i = zeros(N, 1);
    e_i(idfault) = 1;

    % solving for If with the equation on p. 22
    Eeq = e_i' * linsolve(Y,Iint);
    Zeq = e_i' * linsolve(Y,e_i);
    If = Eeq ./ (Zeq + Zf);

    % solving for Vf
    Voc = linsolve(Y, Iint);
    Vf = Voc - If * linsolve(Y, e_i);
end