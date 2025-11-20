%% Q1
mdl = 'q1cd';
load_system(mdl);
set_param(mdl, 'StopTime', '30');

%% Run simulation
out = sim(mdl);

%% q1c

plot_scope(out, 'deltaF', 'Time [s]', '\Delta f [Hz]', 'Part 1: Frequency deviation', 'part1c_freq.png');

%% q1d - Powers

plot_multiple_scopes(out, {'Pm1','Pm2'}, 'Time [s]', '\Delta P_m [MW]', 'Part 1d: Mechanical powers', 'part1d_powers.png', true);

%% Q1e
mdl = 'q1e';
load_system(mdl);
set_param(mdl, 'StopTime', '120');

%% Run simulation
out = sim(mdl);

%% q1e - Frequency

plot_scope(out, 'deltaF', 'Time [s]', '\Delta f [Hz]', 'Part 1(e): Frequency with 1 MW/s limit on Gen 1', 'part1e_freq.png');

%% q1e - Powers

plot_multiple_scopes(out, {'Pm1','Pm2'}, 'Time [s]', '\Delta P_m [MW]', 'Part 1(e): Mechanical powers', 'part1e_powers.png', true);

%% Q2
mdl = 'q2';
load_system(mdl);
set_param(mdl, 'StopTime', '1000');

%% Run simulation
out = sim(mdl);

%% q2 - Frequency

plot_scope(out, 'deltaF', 'Time [s]', '\Delta f [Hz]', 'Part 2: Frequency with AGC', 'part2_freq.png');

%% q2 - Powers

plot_multiple_scopes(out, {'Pm1','Pm2','POverall'}, 'Time [s]', '\Delta P_m [MW]', 'Part 2: Mechanical powers', 'part2_powers.png', true);

%% Q3
mdl = 'ass3q3';
load_system(mdl);
set_param(mdl, 'StopTime', '100');

%% Run simulation
out = sim(mdl);

%% q3 - Frequency

plot_scope(out, 'deltaF', 'Time [s]', '\Delta f [Hz]', 'Part 3: Frequency with lower M and D', 'part3_freq.png');

%% q3 - Powers

plot_multiple_scopes(out, {'Pm1','Pm2','POverall'}, 'Time [s]', '\Delta P_m [MW]', 'Part 3: Mechanical powers', 'part3_powers.png', true);

%% Q4 - No control
mdl = 'ass3q4';
load_system(mdl);
set_param(mdl, 'StopTime', '1200');

%% Run simulation
out = sim(mdl);

%% q4 - Frequency

plot_scopes_separate(out, {'deltaF', 'deltaF1'}, 'Time [s]', '\Delta f [Hz]', 'Part 4: Frequency deviation', 'part4_freq');

%% q4 - System 1 Powers

plot_multiple_scopes(out, {'Pm1','Pm2','POverall'}, 'Time [s]', '\Delta P_m [MW]', 'Part 4: System 1 Mechanical powers', 'part4_sys1_powers.png', true);

%% q4 - System 2 Powers

plot_multiple_scopes(out, {'Pm11','Pm21','POverall1'}, 'Time [s]', '\Delta P_m [MW]', 'Part 4: System 2 Mechanical powers', 'part4_sys2_powers.png', true);

%% q4 - Tie Line Power

plot_scope(out, 'PTieLine', 'Time [s]', 'P_{tie} [MW]', 'Part 4: Tie line power', 'part4_tieline.png');

%% Q4 - With control
mdl = 'ass3q4Control';
load_system(mdl);
set_param(mdl, 'StopTime', '1200');

%% Run simulation
out = sim(mdl);

%% q4 Control - Frequency

plot_scopes_separate(out, {'deltaF', 'deltaF1'}, 'Time [s]', '\Delta f [Hz]', 'Part 4 (Control): Frequency deviation', 'part4c_freq');

%% q4 Control - System 1 Powers

plot_multiple_scopes(out, {'Pm1','Pm2','POverall'}, 'Time [s]', '\Delta P_m [MW]', 'Part 4 (Control): System 1 Mechanical powers', 'part4c_sys1_powers.png', true);

%% q4 Control - System 2 Powers

plot_multiple_scopes(out, {'Pm11','Pm21','POverall1'}, 'Time [s]', '\Delta P_m [MW]', 'Part 4 (Control): System 2 Mechanical powers', 'part4c_sys2_powers.png', true);

%% q4 Control - Tie Line Power

plot_scope(out, 'PTieLine', 'Time [s]', 'P_{tie} [MW]', 'Part 4 (Control): Tie line power', 'part4c_tieline.png');
