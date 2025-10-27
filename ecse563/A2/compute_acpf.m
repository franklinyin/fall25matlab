function [Pf_MW, Sf_MVA] = compute_acpf(Vm, delta, Y, nfrom, nto, Sbase)
    % Compute AC power flows for each transmission line
    % Inputs:
    %   Vm - voltage magnitudes (p.u.)
    %   delta - voltage angles (radians)
    %   Y - admittance matrix
    %   nfrom, nto - line connection vectors
    %   Sbase - base power (MVA)
    % Outputs:
    %   Pf_MW - active power flows (MW)
    %   Sf_MVA - apparent power flows (MVA)
    
    % Construct complex voltages
    Vcomplex = Vm .* exp(1j * delta);
    
    % Get power flow per line
    nl = length(nfrom);
    Pf_MW = zeros(nl,1);
    Sf_MVA = zeros(nl,1);

    for k = 1:nl
        i = nfrom(k); 
        j = nto(k);
        Vi = Vcomplex(i); 
        Vj = Vcomplex(j);

        y_series = -Y(i,j);
        Iij = (Vi - Vj) * y_series;   
        Sij_pu = Vi * conj(Iij);

        Pf_MW(k) = real(Sij_pu) * Sbase;
        Sf_MVA(k) = abs(Sij_pu) * Sbase;
    end
end