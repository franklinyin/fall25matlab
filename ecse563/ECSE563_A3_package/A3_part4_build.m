
function mdl = A3_part4_build()
% Two-area LFC with optional ACE-based secondary control
evalin('base','A3_params');
mdl = 'A3_part4_build';
new_system(mdl); open_system(mdl);
set_param(mdl,'StopTime','1500');

% ---- Helper to build one "Area" as a subsystem ----
function add_area(areaName, xpos)
    add_block('simulink/Ports & Subsystems/Subsystem', [mdl '/' areaName], 'Position', [xpos 80 xpos+420 420]);
    open_system([mdl '/' areaName]);
    % Inside area: same as Part 2 single-area but with inports/outports
    add_block('simulink/Signal Routing/In1', [mdl '/' areaName '/dP_load'], 'Position', [25 210 55 230]); % MW step (external)
    add_block('simulink/Signal Routing/In1', [mdl '/' areaName '/Pref_bus'], 'Position', [25 110 55 130]); % AGC MW bias in
    add_block('simulink/Signal Routing/Out1', [mdl '/' areaName '/f_Hz'], 'Position', [370 140 400 160]);
    add_block('simulink/Signal Routing/In1', [mdl '/' areaName '/Ptie_in'], 'Position', [25 260 55 280]); % +MW leaving area
    
    % Internals (same as single-area)
    add_block('simulink/Math Operations/Gain', [mdl '/' areaName '/Droop1'], 'Gain', '-invR1', 'Position', [130 60 190 90]);
    add_block('simulink/Math Operations/Gain', [mdl '/' areaName '/Droop2'], 'Gain', '-invR2', 'Position', [130 290 190 320]);
    add_block('simulink/Continuous/Transfer Fcn', [mdl '/' areaName '/Gov1'], 'Numerator','[1]','Denominator','[TG1 1]','Position',[220 60 280 90]);
    add_block('simulink/Continuous/Transfer Fcn', [mdl '/' areaName '/Turb1'], 'Numerator','[1]','Denominator','[TCH1 1]','Position',[310 60 370 90]);
    add_block('simulink/Continuous/Transfer Fcn', [mdl '/' areaName '/Gov2'], 'Numerator','[1]','Denominator','[TG2 1]','Position',[220 290 280 320]);
    add_block('simulink/Continuous/Transfer Fcn', [mdl '/' areaName '/Turb2'], 'Numerator','[1]','Denominator','[TCH2 1]','Position',[310 290 370 320]);
    add_block('simulink/Math Operations/Sum', [mdl '/' areaName '/Sum1'], 'Inputs','++','Position',[80 60 110 90]);
    add_block('simulink/Math Operations/Sum', [mdl '/' areaName '/Sum2'], 'Inputs','++','Position',[80 290 110 320]);
    add_block('simulink/Math Operations/Sum', [mdl '/' areaName '/SumP'], 'Inputs','++-+-','Position',[220 170 280 210]);
    add_block('simulink/Math Operations/Gain', [mdl '/' areaName '/D_gain'], 'Gain', 'D', 'Position', [220 230 280 260]);
    add_block('simulink/Continuous/Gain', [mdl '/' areaName '/1_over_M'], 'Gain', '1/M', 'Position', [310 170 370 200]);
    add_block('simulink/Continuous/Integrator', [mdl '/' areaName '/Int_f'], 'Position',[310 210 370 240]);
    add_block('simulink/Sinks/Scope', [mdl '/' areaName '/Scope_f'], 'Position',[310 250 370 280]);
    add_line([mdl '/' areaName],'Int_f/1','Scope_f/1');

    % Pref split handled at top level; Pref_bus enters here and goes equally to both sums
    add_block('simulink/Signal Routing/From', [mdl '/' areaName '/Pref1'], 'GotoTag', areaName, 'Position', [60 110 80 130]);
    add_block('simulink/Signal Routing/From', [mdl '/' areaName '/Pref2'], 'GotoTag', areaName, 'Position', [60 340 80 360]);

    % Wire area internals
    add_line([mdl '/' areaName], 'Pref1/1', 'Sum1/1'); add_line([mdl '/' areaName], 'Droop1/1','Sum1/2');
    add_line([mdl '/' areaName], 'Pref2/1', 'Sum2/1'); add_line([mdl '/' areaName], 'Droop2/1','Sum2/2');
    add_line([mdl '/' areaName], 'Sum1/1','Gov1/1'); add_line([mdl '/' areaName], 'Gov1/1','Turb1/1');
    add_line([mdl '/' areaName], 'Sum2/1','Gov2/1'); add_line([mdl '/' areaName], 'Gov2/1','Turb2/1');
    % Pm sum -> swing
    add_line([mdl '/' areaName], 'Turb1/1','SumP/1'); add_line([mdl '/' areaName], 'Turb2/1','SumP/2');
    add_line([mdl '/' areaName], 'dP_load/1','SumP/3'); add_line([mdl '/' areaName], 'D_gain/1','SumP/4');
    add_line([mdl '/' areaName], 'Ptie_in/1','SumP/5');
    add_line([mdl '/' areaName], 'SumP/1','1_over_M/1'); add_line([mdl '/' areaName], '1_over_M/1','Int_f/1');
    add_line([mdl '/' areaName], 'Int_f/1','Droop1/1'); add_line([mdl '/' areaName], 'Int_f/1','Droop2/1');
    add_line([mdl '/' areaName], 'Int_f/1','f_Hz/1');
    close_system([mdl '/' areaName]);
end

% Build Area1 and Area2 subsystems
add_area('Area1', 50);
add_area('Area2', 650);

% Top-level AGC blocks (frequency-only integrators first)
add_block('simulink/Math Operations/Gain', [mdl '/Ki1'], 'Gain','Ki_AGC','Position',[120 20 180 40]);
add_block('simulink/Math Operations/Gain', [mdl '/Ki2'], 'Gain','Ki_AGC','Position',[720 20 780 40]);
add_block('simulink/Math Operations/Gain', [mdl '/Neg1'], 'Gain','-1','Position',[50 20 110 40]);
add_block('simulink/Math Operations/Gain', [mdl '/Neg2'], 'Gain','-1','Position',[650 20 710 40]);
add_block('simulink/Continuous/Integrator', [mdl '/Int_u1'], 'Position',[190 20 250 40]);
add_block('simulink/Continuous/Integrator', [mdl '/Int_u2'], 'Position',[790 20 850 40]);
add_block('simulink/Math Operations/Gain', [mdl '/a1_A1'], 'Gain','a1','Position',[260 20 320 40]);
add_block('simulink/Math Operations/Gain', [mdl '/a2_A1'], 'Gain','a2','Position',[330 20 390 40]);
add_block('simulink/Math Operations/Gain', [mdl '/a1_A2'], 'Gain','a1','Position',[860 20 920 40]);
add_block('simulink/Math Operations/Gain', [mdl '/a2_A2'], 'Gain','a2','Position',[930 20 990 40]);
add_block('simulink/Signal Routing/Goto', [mdl '/BusA1'], 'GotoTag','Area1','Position',[400 20 460 40]);
add_block('simulink/Signal Routing/Goto', [mdl '/BusA2'], 'GotoTag','Area2','Position',[1000 20 1060 40]);

% Wire AGC (frequency only)
add_line(mdl,'Area1/f_Hz','Neg1/1'); add_line(mdl,'Neg1/1','Ki1/1'); add_line(mdl,'Ki1/1','Int_u1/1');
add_line(mdl,'Int_u1/1','a1_A1/1'); add_line(mdl,'Int_u1/1','a2_A1/1');
add_line(mdl,'a1_A1/1','BusA1/1'); add_line(mdl,'a2_A1/1','BusA1/1');

add_line(mdl,'Area2/f_Hz','Neg2/1'); add_line(mdl,'Neg2/1','Ki2/1'); add_line(mdl,'Ki2/1','Int_u2/1');
add_line(mdl,'Int_u2/1','a1_A2/1'); add_line(mdl,'Int_u2/1','a2_A2/1');
add_line(mdl,'a1_A2/1','BusA2/1'); add_line(mdl,'a2_A2/1','BusA2/1');

% Load steps
add_block('simulink/Sources/Step', [mdl '/Load1'], 'Time','0','Before','0','After','DeltaP','Position',[50 440 80 460]);
add_block('simulink/Sources/Step', [mdl '/Load2'], 'Time','0','Before','0','After','0','Position',[650 440 680 460]);
add_line(mdl,'Load1/1','Area1/dP_load');
add_line(mdl,'Load2/1','Area2/dP_load');

% Tie-line: d(delta)/dt = 2*pi*(f1 - f2); P12 = (Psys/Xline) * delta
add_block('simulink/Math Operations/Sum', [mdl '/fDiff'], 'Inputs','+-','Position',[530 170 560 190]);
add_block('simulink/Math Operations/Gain', [mdl '/TwoPi'], 'Gain', '2*pi', 'Position',[570 170 630 190]);
add_block('simulink/Continuous/Integrator', [mdl '/Int_delta'], 'InitialCondition','0','Position',[640 170 700 190]);
add_block('simulink/Math Operations/Gain', [mdl '/T12'], 'Gain', 'Psys/0.1', 'Position',[710 170 770 190]); % 0.1 pu line reactance
add_block('simulink/Math Operations/Gain', [mdl '/Neg'], 'Gain','-1','Position',[780 170 820 190]);

% Connect tie-line signals
add_line(mdl,'Area1/f_Hz','fDiff/1'); add_line(mdl,'Area2/f_Hz','fDiff/2');
add_line(mdl,'fDiff/1','TwoPi/1'); add_line(mdl,'TwoPi/1','Int_delta/1');
add_line(mdl,'Int_delta/1','T12/1');
add_line(mdl,'T12/1','Area1/Ptie_in');  % +MW leaving Area1
add_line(mdl,'T12/1','Neg/1'); add_line(mdl,'Neg/1','Area2/Ptie_in'); % -MW entering Area2

% Scopes
add_block('simulink/Sinks/Scope', [mdl '/Scope_Ptie'], 'Position',[840 170 880 200]);
add_line(mdl,'T12/1','Scope_Ptie/1');

save_system(mdl);
end
