
%% Part 4 — Two-area tests
A3_params;
mdl4 = A3_part4_build();
set_param(mdl4,'StopTime','1500');

% Case A: frequency-only integral control (NO ACE) — tie-line flow will not return to schedule
out4A = sim(mdl4);

% Enable ACE: ACE1 = B*f1 + P12; ACE2 = B*f2 - P12
% Build minimal ACE logic at the top level by replacing Ki inputs
open_system(mdl4);
% Add bias B * f blocks and sum with tie-line for each area
add_block('simulink/Math Operations/Gain', [mdl4 '/B1'], 'Gain','(D+invR1+invR2)','Position',[50 500 110 520]);
add_block('simulink/Math Operations/Gain', [mdl4 '/B2'], 'Gain','(D+invR1+invR2)','Position',[650 500 710 520]);
add_block('simulink/Math Operations/Sum', [mdl4 '/ACE1sum'], 'Inputs','++','Position',[120 500 150 520]);
add_block('simulink/Math Operations/Sum', [mdl4 '/ACE2sum'], 'Inputs','+-','Position',[720 500 750 520]);

add_line(mdl4,'Area1/f_Hz','B1/1'); add_line(mdl4,'B1/1','ACE1sum/1'); add_line(mdl4,'T12/1','ACE1sum/2');
add_line(mdl4,'Area2/f_Hz','B2/1'); add_line(mdl4,'B2/1','ACE2sum/1'); add_line(mdl4,'T12/1','ACE2sum/2');

% Rewire Ki blocks to use -ACE instead of -f
delete_line(mdl4,'Area1/f_Hz/1','Neg1/1'); delete_line(mdl4,'Area2/f_Hz/1','Neg2/1');
add_block('simulink/Math Operations/Gain', [mdl4 '/NegACE1'], 'Gain','-1','Position',[160 500 180 520]);
add_block('simulink/Math Operations/Gain', [mdl4 '/NegACE2'], 'Gain','-1','Position',[760 500 780 520]);
add_line(mdl4,'ACE1sum/1','NegACE1/1'); add_line(mdl4,'NegACE1/1','Ki1/1');
add_line(mdl4,'ACE2sum/1','NegACE2/1'); add_line(mdl4,'NegACE2/1','Ki2/1');

save_system(mdl4);
out4B = sim(mdl4);
disp('Part 4 simulations complete: out4A (no ACE), out4B (with ACE).');
