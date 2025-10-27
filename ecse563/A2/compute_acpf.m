% Q5 helper function
function Sij_complex = compute_acpf(Vm, delta, Y, nfrom, nto, Sbase)
    % Helper function for Q5 to compute AC power flows for each transmission line
    % Inputs:
    %   Vm - voltage magnitudes (p.u.)
    %   delta - voltage angles (radians)
    %   Y - admittance matrix
    %   nfrom, nto - line connection vectors
    %   Sbase - base power (MVA)
    % Output:
    %   Sij_complex - complex power flows (MVA)
    
    % Construct complex voltages
    Vcomplex = Vm .* exp(1j * delta);
    
    % Get power flow per line
    nl = length(nfrom);
    Sij_complex = zeros(nl,1);

    for k = 1:nl
        i = nfrom(k); 
        j = nto(k);
        Vi = Vcomplex(i); 
        Vj = Vcomplex(j);

        y_series = -Y(i,j);
        Iij = (Vi - Vj) * y_series;   
        Sij_pu = Vi * conj(Iij);

        Sij_complex(k) = Sij_pu * Sbase;
    end
end