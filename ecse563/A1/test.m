% Load IEEE 9-bus system data
run('ieee9_A1.m');


%% 1. 
disp('Part 1:');
disp('Testing admittance function...');
Y = admittance(nfrom, nto, r, x, b);

disp('Admittance matrix (Y):');
disp(Y);

% V = IZ 
V = linsolve(Y, Iint);

disp('Voltage (V):');
disp(V);

magnitude = abs(V);
angle_rad = angle(V);
angle_deg = angle_rad * (180/pi);

polar_coordinates_rad = [magnitude, angle_rad];
disp('  Magnitude   Angle (in radian)');
disp(polar_coordinates_rad);

polar_coordinates_deg = [magnitude, angle_deg];
disp('  Magnitude   Angle (in degree)');
disp(polar_coordinates_deg);

%% 2.
disp('Part 2:');
Z = impedance(nfrom, nto, r, x, b);  


disp('Impedance matrix (Z):');
disp(Z);

%% 3. 
disp('Part 3:');
disp('Testing on Node 3 for example');
idf = 3;
Zf = 0;

[If, Vf] = fault(Y, Iint, idf, Zf); 

disp('Fault current at node 3 (If):');
disp(If);
disp('Node voltages during fault (Vf):');
disp(Vf);

% then we compute for all nodes (lines), and plot them
figure;
hold on;
colors = lines(9);

for idf = 1:9
    [If, Vf] = fault(Y, Iint, idf, Zf);
    Vabs = abs(Vf);
    plot(Vabs, 'Color', colors(idf, :), 'DisplayName', sprintf('idf = %d, If = %.4f + %.4fi', idf, real(If), imag(If)));
    text(length(Vabs), Vabs(end), sprintf('idf = %d, If = %.4f + %.4fi', idf, real(If), imag(If)), 'VerticalAlignment', 'bottom');
end

legend show;
axis padded; 
xlim([0, length(Vabs) + 5]);
legend('Location', 'eastoutside');
xlabel('Node');
ylabel('|V|');
title('Q3');



%% 4. Test generalized Thevenin equivalent function
disp('Part 4:');
disp('a)');
id_thev = [3, 5];  

[Eeq, Zeq] = genthevenin(Y, Iint, id_thev);


disp(['For nodes ', num2str(id_thev)]);
disp('Thevenin equivalent voltages (Eeq):');
disp(Eeq);
disp('Thevenin equivalent impedance matrix (Zeq):');
disp(Zeq);



disp('b)');
id_thev = [9, 4]; 

[Eeq, Zeq] = genthevenin(Y, Iint, id_thev);

% Display the Thevenin equivalent voltages and impedance matrix
disp(['For nodes ', num2str(id_thev)]);
disp('Thevenin equivalent voltages (Eeq):');
disp(Eeq);
disp('Thevenin equivalent impedance matrix (Zeq):');
disp(Zeq);

% some additional variables for the plot:
disp('I:');
disp(Iint(9)); % the 9th channel connecting 9-4

disp('delta y:');
disp(Y(9,4)); % Y(4,9) works too


%% ==========part 5============
%% 5. Test generalized fault calculation function
disp('Part 5:');
disp('a)');
% a) 
disp('Testing generalized fault calculation function...');

YN = Y;
% Create fault admittance matrix by removing line between nodes 8 and 9
YF = outage_admittance_helper(YN, nfrom, nto, r, x, b, 8, 9);

disp('YF:');
disp(YF);


[If, Vf] = fault(YF, Iint, 8, 0);


IintF = Iint;

idN = [8, 9]';
idF = [8, 9]';

[IT, VNF] = genfault(YN, YF, Iint, IintF, idN, idF)

print_magnitudes_angles(IT, VNF, idN, idF);

% b) 
disp('b)');

idN = [1];
idF = [5];

[IT, VNF] = genfault(YN, YN, Iint, IintF, idN, idF)

print_magnitudes_angles(IT, VNF, idN, idF);

% c) 
disp('c)');
idN = [3 5]';
idF = [7 4]'; 

[IT, VNF] = genfault(YN, YN, Iint, IintF, idN, idF)

print_magnitudes_angles(IT, VNF, idN, idF);



% d) 
disp('d)');



% Load IEEE 24-bus system data
run('ieee24_A1.m');

% Create two identical IEEE 24-bus systems
Y = admittance(nfrom, nto, r, x, b);
YN_24 = Y;  % First system (healthy network)
IintN_24 = Iint;  % Current sources for first system

YF_24 = Y;  % Second system (fault network) 
IintF_24 = Iint;  % Current sources for second system

idN = [7 13 23]';
idF = [3 15 17]'; 

[IT_24, VNF_24] = genfault(YN_24, YN_24, IintN_24, IintN_24, idN, idF)

print_magnitudes_angles(IT_24, VNF_24, idN, idF);



disp('All tests completed successfully.');
