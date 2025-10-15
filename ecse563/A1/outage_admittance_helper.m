function YF = outage_admittance_helper(YN, nfrom, nto, r, x, b, outage_from, outage_to)
    % creates fault admittance matrix by removing a line
    %
    % Inputs:
    %   YN - original healthy network admittance matrix
    %   nfrom - vector of sending buses
    %   nto - vector of receiving buses
    %   r, x, b - vectors of resistance, reactance and shunt susceptance
    %   outage_from - from which node of the line to remove
    %   outage_to - to which node of the line to remove
    %
    % Output:
    %   YF - Fault admittance matrix with specified line removed
    
    % find the line index for the specified outage
    line_idx = find((nfrom == outage_from & nto == outage_to) | (nfrom == outage_to & nto == outage_from));
    
    % Initialize fault admittance matrix as copy of healthy network
    YF = YN;
    
    % calculate line admittance and shunt susceptance
    y_line = 1/(r(line_idx) + 1i*x(line_idx));  % line admittance
    b_shunt = 1i*b(line_idx)/2;  % half of shunt susceptance
    
    % Remove the line admittance from the matrix
    YF(outage_from, outage_to) = YF(outage_from, outage_to) + y_line;  % add back the negative admittance (remove it)
    YF(outage_to, outage_from) = YF(outage_to, outage_from) + y_line;  % symmetry
    YF(outage_from, outage_from) = YF(outage_from, outage_from) - y_line - b_shunt;  % remove from diagonal
    YF(outage_to, outage_to) = YF(outage_to, outage_to) - y_line - b_shunt;  % remove from diagonal

end
