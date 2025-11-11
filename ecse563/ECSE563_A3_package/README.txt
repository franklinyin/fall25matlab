
ECSE 563 — Assignment 3 (Generation & Frequency Control)

FILES
-----
A3_params.m        : Base parameters (units: MW, Hz, s)
A3_part1_build.m   : Builds Part 1 single-area model (primary only). Rate Limiter block included but commented.
A3_part2_build.m   : Builds Part 2 model (adds AGC integral + participation).
A3_part3_build.m   : Same as Part 2; run with M=0.6*M, D=0.5*D for Part 3.
A3_part4_build.m   : Two-area model (frequency-only integral by default).
A3_run_all.m       : Script to build & run Parts 1–3 and save figures.
run_part4.m        : Script to run Part 4 (without ACE, then with ACE).

REPORT
------
report.tex         : LaTeX report with filled-in calculations; compile after running A3_run_all to embed figures.

HOW TO RUN
----------
1) Open MATLAB, change folder to this directory.
2) Run:  A3_params; A3_run_all
3) For Part 4: run:  A3_params; run_part4

TIPS
----
- All signals are in MW and Hz. The swing equation uses M in MW*s/Hz.
- If scopes are empty, enable data logging or open the scopes while simulating.
- If Part 4 oscillates, reduce Ki_AGC or add small lag/washout in ACE.
