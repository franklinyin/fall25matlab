function print_magnitudes_angles(IT, VNF, idN, idF)
    % helper function printing magnitudes and angles of IT and VNF vectors
    %
    % Inputs:
    %   IT - Tie-line currents vector
    %   VNF - Healthy network node voltages with fault vector
    %   idN - Healthy node indices
    %   idF - Fault node indices
    
    % print IT magnitudes and angles
    fprintf('IT:\n');
    for i = 1:length(IT)
        fprintf('line N%d-F%d: |I| = %.4f, angle = %8.2f°\n', idN(i), idF(i), abs(IT(i)), angle(IT(i)) * 180/pi);
    end
    
    % Print VNF magnitudes and angles
    fprintf('\nVNF:\n');
    for i = 1:length(VNF)
        fprintf('Line %2d: |V| = %.4f, angle = %8.2f°\n', i, abs(VNF(i)), angle(VNF(i)) * 180/pi);
    end
    
    fprintf('\n');
end
