
# ECSE 563 A3 — Simulink Build Cheatsheet (Click-by-click)

This is **painfully specific** so you can rebuild everything without knowing Simulink.

> **Notation**: In Simulink R2021b+ the main toolbar is called the *Toolstrip*. If your layout differs, use `Home ▸ Library Browser` to find blocks.

---

## A. Open Simulink and set up the folder
1. Launch **MATLAB**.
2. In the **Current Folder** panel (left), **Right-click** and choose **New ▸ Folder**, name it `ECSE563_A3`.
3. **Drag** the downloaded files into this folder: `A3_params.m`, `A3_part1_build.m`, `A3_part2_build.m`, `A3_part3_build.m`, `A3_part4_build.m`, `A3_run_all.m`.
4. In the **Command Window**, type:
   ```matlab
   cd ECSE563_A3
   A3_params
   ```

---

## B. Part 1 model (primary control only)
You can either run the builder script or place blocks by hand.

### Option 1 — Build by script (recommended)
1. In the Command Window: `mdl1 = A3_part1_build;`
2. A new diagram named **A3\_part1\_build** opens with all blocks laid out.
3. On the top **Simulation** tab, set **Stop Time** to `30` seconds.
4. Click **Run** (green triangle). Open **Scope\_f** to see frequency and **Scope\_Pm** to see Gen1/Gen2 mechanical powers.

### Option 2 — Build by hand
1. **Home ▸ New ▸ Model** → a blank canvas opens.
2. **Library Browser** (left side) → expand **Simulink**:
   - **Sources** → drag **Step** block; double-click it and set *After* = `100` (MW), *Time* = `0`.
   - **Math Operations** → drag two **Gain** blocks; rename them `Droop1` and `Droop2`. Double-click and set **Gain** to `-166.6667` and `-83.3333` (these are `-1/R1` and `-1/R2` in MW/Hz).
   - **Continuous** → drag two **Transfer Fcn** blocks for each generator: name them `Gov1`,`Turb1` and `Gov2`,`Turb2`. Set **Denominator** to `[TG 1]` (`[0.25 1]`, `[3 1]`, `[0.5 1]`, `[8 1]`).
   - **Math Operations** → drag two **Sum** blocks, set **Inputs** to `++` (click the block, in the property dialog set `List of signs` to `++`). These sum Pref and Droop for each unit.
   - **Math Operations** → drag a **Sum** block named `SumP`, set **Inputs** to `++-+` (Pm1 + Pm2 − Step + D*f).
   - **Math Operations** → drag a **Gain** block named `D_gain`, set **Gain** = `20` (MW/Hz).
   - **Continuous** → drag **Gain** (`1_over_M`, with **Gain**=`1/700`) and an **Integrator** (`Int_f`) for the swing equation.
   - **Sinks** → drag two **Scope** blocks (`Scope_f`, `Scope_Pm`), and two **To Workspace** blocks (names: `f_out` with *Variable name* `f_Hz`, `Pm_out` with `Pm_MW`).
   - **Signal Routing** → drag a **Mux** (2 inputs) to combine `Pm1` and `Pm2` for plotting.
3. Wiring (use click-drag from output arrow to input triangle):
   - `Int_f` → `Droop1` and `Droop2` (this is *f*).
   - `Droop1` + `Pref1(=0)` → `Sum1` → `Gov1` → `Turb1` → (label this line `Pm1`).
   - `Droop2` + `Pref2(=0)` → `Sum2` → `Gov2` → `Turb2` → (label this line `Pm2`).
   - `Pm1` and `Pm2` into `SumP` (+ +); `Step` into `SumP` (−); `D_gain` into `SumP` (+).
   - `Int_f` ← `1_over_M` ← `SumP`.
   - `Int_f` → `Scope_f` and `f_out`; `Mux` takes `Pm1` & `Pm2` → `Scope_Pm` and `Pm_out`.
4. On the **Simulation** tab set **Stop Time = 30** and click **Run**.

**Steady-state check**: From **Scope\_f**, the final value should be about **−0.370 Hz**, i.e. **59.63 Hz**.

---

## C. Part 1(e): add a ramp-rate limit to Gen 1
1. **Library Browser ▸ Discontinuities** → drag **Rate Limiter**. Drop it between `Turb1` and `Pm1`.
2. Double-click **Rate Limiter**, set **Rising slew limit = 1** and **Falling slew limit = 1** (MW/s). Sample time `0` (continuous).
3. Increase **Stop Time** to `120` s. **Run**. The frequency nadir will deepen (≈ −0.75 Hz around 30 s).

*(If you used the script build: un-comment the block via right-click → **Comment Out** to toggle.)*

---

## D. Part 2: add secondary (integral) control with participation
1. Use the script: `mdl2 = A3_part2_build;` (preferred), or by hand:
   - **Math Operations** → place **Gain** `Ki` (= `1.0`), **Gain** `-1` then **Integrator** `Int_u` (this integrates `-f` in MW).
   - Split `Int_u` with **Gain** `a1=0.75` to **Pref1**, and **Gain** `a2=0.25` to **Pref2`.
   - Wire `Pref1` to `Sum1` and `Pref2` to `Sum2` (replacing the zeros from Part 1).
2. Set **Stop Time = 900** s. **Run**. Frequency should recover near 0 Hz by 10–15 minutes; the final power sharing is **75 MW** (Gen1) and **25 MW** (Gen2).

*Tip:* If you see oscillations, try reducing `Ki` to `0.8` or add a small lead (e.g., replace `Ki` with a first-order lag `Ki/(1+3s)`).

---

## E. Part 3: lower inertia and load damping
1. In the MATLAB Command Window, before running the model, type:
   ```matlab
   M = 0.6*M; D = 0.5*D;
   ```
2. Use `A3_part3_build;` and **Run** for 900 s.
3. Expect a deeper nadir and slightly more oscillatory behaviour with the same `Ki` — discuss whether `Ki` should be re-tuned or fast frequency response added.

---

## F. Part 4: two-area interconnection (with and without ACE)
1. Build via script: `mdl4 = A3_part4_build;`
   - The diagram contains **Area1** and **Area2** subsystems (copies of Part 2), a tie-line integrator (Δδ), and a gain `Psys/0.1` to implement a line with `X = 0.1 p.u.`.
   - `Load1` applies the 100 MW step to **Area1** only.
2. **Run** (1500 s). With **frequency-only** secondary, frequency errors go to zero in both areas **but** the tie-line settles to a nonzero value (Area 1 imports).
3. Now enable ACE (Area Control Error) by executing the lines in `run_part4.m` or following the in-model notes:
   - ACE\_1 = B·f\_1 + P\_12; ACE\_2 = B·f\_2 − P\_12 with `B ≈ 270 MW/Hz`.
   - Rewire the `Ki` blocks to integrate `−ACE` instead of `−f` (the script does this automatically).
4. **Run** again. With ACE, both **frequency** and **tie-line deviation** return to zero. If you see slow oscillations, reduce `Ki` or add a small washout on ACE (first-order high-pass).

---

## G. Export results and plots
- Use the provided **A3\_run\_all.m** to run Parts 1–3 and save figures (`part1_freq.png`, etc.).
- For Part 4, log `Int_delta` (angle) and `T12` (power) plus both area frequencies on a scope and use **Save Data to Workspace** for plotting.

---

## H. Quick sanity numbers (so you can answer viva):
- $R_1 = 0.006$ Hz/MW, $R_2 = 0.012$ Hz/MW; on 2000 MW base: $12$ and $24$ Hz/p.u. → $0.20$ and $0.40$ p.u./p.u.
- $M = 700$ MW·s/Hz ($H = 10.5$ s on 2000 MW base).
- Primary-only steady frequency: $59.63$ Hz (error $pprox -0.370$ Hz).
- Primary sharing at steady state: Gen1 $pprox 61.7$ MW, Gen2 $pprox 30.9$ MW, load relief $pprox -7.4$ MW.
- With AGC ($K_I pprox 1.0$ MW/(Hz·s)), frequency $ightarrow 60$ Hz in 10–15 min; final sharing 75/25 MW as per participation.
- With lower $M$ and $D$: deeper nadir, ROCOF faster, more oscillation — consider retuning $K_I$ or adding synthetic inertia.
- Two-area, no ACE: frequency ok, tie-line not ok (Area1 imports $\sim$50 MW for symmetric areas). With ACE: both ok.

Good luck!
