
%% ECSE563 - Assignment 3 - Part 1(c)
% Programmatic build of the primary frequency control model (Bouffard notation, rad/s).
% Run this script in MATLAB. It will create/open the Simulink model and set all parameters.
%
% Model contents:
%   - Two generator control paths: droop -> governor (1/(1+s*TG)) -> turbine (1/(1+s*TCH))
%   - Plant (generator + load): Δω(s) = (ΔPm1+ΔPm2 - ΔPL - D*Δω) / (M*s + D)
%   - Disturbance: Step ΔPL = +100 MW at t=0 (subtracting at plant input)
%   - Scopes for Δω (rad/s & Hz) and ΔPm1, ΔPm2
%
% Sources:
%   Bouffard slides (ECSE563_frequency-control_2025.pdf): swing eqn / M=ω0 J; droop R = Δω/ΔPe;
%   Assignment details (ECSE563_assignment_3_F2025.pdf): TG/TCH, ratings, Wk0, D in MW/Hz.
%
% After build, click 'Run' in Simulink or call: sim('ecse563_part1c_model')

%% Parameters (numerical values from 1(a) and 1(b))
f0   = 60;                        % Hz
w0   = 2*pi*f0;                   % rad/s
Wk0J = 21e9;                      % J  (stored kinetic energy)
M    = (2*Wk0J)/w0/1e6;           % MW*s^2/rad  (convert J->MW*s by /1e6 implicitly via units)
M    = 111.4;                     % lock to computed value to avoid roundoff

D_Hz = 20;                        % MW/Hz (load frequency sensitivity)
D    = D_Hz/(2*pi);               % MW/(rad/s)

% Droop constants from 1(a)
R1_HzMW = 0.006;                  % Hz/MW
R2_HzMW = 0.012;                  % Hz/MW
R1w     = 2*pi*R1_HzMW;           % rad/s per MW
R2w     = 2*pi*R2_HzMW;           % rad/s per MW

% Governor / turbine time constants
TG1  = 0.25;   TCH1 = 3.0;
TG2  = 0.50;   TCH2 = 8.0;

% Disturbance (load step)
dPL  = 100;                       % MW, positive for load increase (subtract at plant input)

%% Create model
mdl = 'ecse563_part1c_model';
if bdIsLoaded(mdl); close_system(mdl,0); end
new_system(mdl); open_system(mdl);

x = 30;  y = 30; dx = 140; dy = 70;  % layout helpers

%% Blocks: Plant
add_block('simulink/Math Operations/Sum', [mdl '/Sum_Plant'], ...
    'Position', [x+400 y+50 x+420 y+130], 'Inputs', '++-');
% Inputs to Sum_Plant: +ΔPm_total, +(-D*Δω), -ΔPL

add_block('simulink/Continuous/Transfer Fcn', [mdl '/Plant_TF'], ...
    'Position', [x+480 y+50 x+600 y+130], 'Numerator', '1', ...
    'Denominator', sprintf('[%g %g]', M, D));
% Δω(s) = (ΔPm_total - ΔPL - D*Δω)/ (M*s + D)

% Feedback path for D*Δω
add_block('simulink/Math Operations/Gain', [mdl '/D_Gain'], ...
    'Position', [x+660 y+160 x+710 y+200], 'Gain', num2str(D));

add_block('simulink/Sinks/Scope', [mdl '/Scope_omega'], ...
    'Position', [x+680 y+40 x+740 y+110]);
set_param([mdl '/Scope_omega'], 'NumInputPorts', '2'); % Δω and Δf

% Δf conversion
add_block('simulink/Math Operations/Gain', [mdl '/ToHz'], ...
    'Position', [x+620 y+100 x+670 y+140], 'Gain', '1/(2*pi)');

%% Blocks: ΔPL Step (subtract at plant sum)
add_block('simulink/Sources/Step', [mdl '/Step_Load'], ...
    'Position', [x+320 y+100 x+350 y+130], 'Time', '0', 'Before', '0', 'After', num2str(dPL));

%% Blocks: Generator 1 path (droop -> governor -> turbine)
add_block('simulink/Math Operations/Gain', [mdl '/Droop1'], ...
    'Position', [x+160 y x+210 y+40], 'Gain', sprintf('-%g', 1/R1w));
% Note: Droop Gain is -1/R (MW per (rad/s))

add_block('simulink/Continuous/Transfer Fcn', [mdl '/Gov1'], ...
    'Position', [x+230 y x+310 y+40], 'Numerator', '1', ...
    'Denominator', sprintf('[%g 1]', TG1));

add_block('simulink/Continuous/Transfer Fcn', [mdl '/Turb1'], ...
    'Position', [x+330 y x+410 y+40], 'Numerator', '1', ...
    'Denominator', sprintf('[%g 1]', TCH1));

%% Blocks: Generator 2 path
add_block('simulink/Math Operations/Gain', [mdl '/Droop2'], ...
    'Position', [x+160 y+120 x+210 y+160], 'Gain', sprintf('-%g', 1/R2w));

add_block('simulink/Continuous/Transfer Fcn', [mdl '/Gov2'], ...
    'Position', [x+230 y+120 x+310 y+160], 'Numerator', '1', ...
    'Denominator', sprintf('[%g 1]', TG2));

add_block('simulink/Continuous/Transfer Fcn', [mdl '/Turb2'], ...
    'Position', [x+330 y+120 x+410 y+160], 'Numerator', '1', ...
    'Denominator', sprintf('[%g 1]', TCH2));

%% Sum ΔPm1 + ΔPm2
add_block('simulink/Math Operations/Sum', [mdl '/Sum_Pm'], ...
    'Position', [x+420 y+60 x+440 y+120], 'Inputs', '++');

%% Scopes for Pm
add_block('simulink/Sinks/Scope', [mdl '/Scope_Pm'], ...
    'Position', [x+420 y+180 x+480 y+240]);
set_param([mdl '/Scope_Pm'], 'NumInputPorts', '2');

%% Lines: wiring
% Δω feedback to droops
add_line(mdl, 'Plant_TF/1', 'ToHz/1');
add_line(mdl, 'ToHz/1', 'Scope_omega/2');
add_line(mdl, 'Plant_TF/1', 'Droop1/1');
add_line(mdl, 'Plant_TF/1', 'Droop2/1');

% Droop1 -> Gov1 -> Turb1 -> Sum_Pm
add_line(mdl, 'Droop1/1', 'Gov1/1');
add_line(mdl, 'Gov1/1', 'Turb1/1');
add_line(mdl, 'Turb1/1', 'Sum_Pm/1');
add_line(mdl, 'Turb1/1', 'Scope_Pm/1');

% Droop2 -> Gov2 -> Turb2 -> Sum_Pm
add_line(mdl, 'Droop2/1', 'Gov2/1');
add_line(mdl, 'Gov2/1', 'Turb2/1');
add_line(mdl, 'Turb2/1', 'Sum_Pm/2');
add_line(mdl, 'Turb2/1', 'Scope_Pm/2');

% Sum_Pm -> Sum_Plant (+)
add_line(mdl, 'Sum_Pm/1', 'Sum_Plant/1');

% Plant_TF output -> Scope_omega (Δω)
add_line(mdl, 'Plant_TF/1', 'Scope_omega/1');

% Δω -> D_Gain -> Sum_Plant (+) second input (as +(-D*Δω))
add_line(mdl, 'Plant_TF/1', 'D_Gain/1');
add_line(mdl, 'D_Gain/1', 'Sum_Plant/2');

% Load step -> Sum_Plant (-) third input
add_line(mdl, 'Step_Load/1', 'Sum_Plant/3');

% Sum_Plant -> Plant_TF
add_line(mdl, 'Sum_Plant/1', 'Plant_TF/1');

%% Clean up & save
set_param(mdl, 'StopTime', '30');
save_system(mdl);
open_system(mdl);
disp('Model ecse563_part1c_model created. Click Run in Simulink or call sim(mdl).');
