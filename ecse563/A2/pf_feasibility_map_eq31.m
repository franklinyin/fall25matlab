function [feas, params] = pf_feasibility_map_eq31(Y, is, busj, Sbase, PvalsMW, QvalsMVAr, V1mag)
% pf_feasibility_map_eq31
% Feasibility of (P,Q) at bus 'busj' using the two-bus solvability inequality (slide eq. (31)).
% Inputs:
%   Y          : full bus admittance matrix (from admittance.m)
%   is         : index of slack bus (sending end of the two-bus)
%   busj       : index of the bus of interest (e.g., 7)
%   Sbase      : MVA base
%   PvalsMW    : vector of P (MW) values to sweep at busj
%   QvalsMVAr  : vector of Q (MVAr) values to sweep at busj
%   V1mag      : sending-end voltage magnitude |V1| (p.u.). Use your spec (e.g., 1.0).
% Outputs:
%   feas       : logical matrix [numel(QvalsMVAr) x numel(PvalsMW)], true if feasible
%   params     : struct with b_eq (p.u. susceptance), y_series, and Yred used

    % ---- 1) Build a series-only Y (remove shunts) so reduction gives the series element cleanly
    Ys = Y;                                   % copy
    Ys(1:size(Y,1)+1:end) = 0;                % zero out the diagonal
    diagSeries = -sum(Ys,2);                  % series-only self-admittance from off-diagonals
    Ys(1:size(Y,1)+1:end) = diagSeries;       % put back on the diagonal

    % ---- 2) Kron reduce {all but [is, busj]} to get a 2x2 network
    keep = [is, busj];
    elim = setdiff(1:size(Ys,1), keep);

    Y11  = Ys(keep, keep);
    Y1e  = Ys(keep, elim);
    Ye1  = Ys(elim, keep);
    Yee  = Ys(elim, elim);

    % Avoid explicit inverse; solve Yee * X = Ye1
    X    = Yee \ Ye1;                         % Yee * X = Ye1
    Yred = Y11 - Y1e * X;                     % Schur complement (Kron reduction)

    % ---- 3) Equivalent series admittance and its susceptance magnitude
    y_series = -Yred(1,2);                    % series element between keep(1)=is and keep(2)=busj
    b_eq     = abs(imag(y_series));           % |Im{y_series}| (p.u.)

    % Safety check
    if b_eq <= eps
        error('Equivalent susceptance is ~0; check network connectivity or data.');
    end

    % ---- 4) Build the feasibility map using ( (2Q/b - V1^2)^2 >= 4/b^2 * (P^2 + Q^2) )
    [PP_MW, QQ_MVAr] = meshgrid(PvalsMW, QvalsMVAr);
    p = PP_MW / Sbase;                        % per unit
    q = QQ_MVAr / Sbase;                      % per unit

    disc = (2*q./b_eq - V1mag.^2).^2 - 4./(b_eq.^2) .* (p.^2 + q.^2);
    feas = disc >= 0;                         % equality is the boundary

    % ---- 5) (Optional) plot with boundary curve
    figure('Color','w'); hold on; box on;
    % background (infeasible) and overlay feasible
    imagesc(PvalsMW, QvalsMVAr, ~feas);       % 1 -> infeasible, 0 -> feasible (blue-ish default)
    colormap([0.8 0.9 1.0; 1.0 0.8 0.9]);     % light blue (infeasible), light pink (feasible)
    set(gca,'YDir','normal');
    % mask to paint feasible region
    h = imagesc(PvalsMW, QvalsMVAr, feas);
    set(h,'AlphaData',0.6);                   % semi-transparent feasible overlay

    % Overlay analytic boundary (disc == 0) as a black curve
    k  = 0.5 * V1mag^2 * b_eq;                % handy constant
    ppu = PvalsMW / Sbase;
    q_boundary_pu = (k^2 - ppu.^2) / (2*k);   % rearranged equality -> q(p)
    plot(PvalsMW, q_boundary_pu * Sbase, 'k-', 'LineWidth', 1.5);

    xlabel('P_7 (MW)'); ylabel('Q_7 (MVAr)');
    title(sprintf('Two-bus solvability map at bus %d  (V_1=%.3f p.u.,  b_{eq}=%.4f p.u.)', ...
            busj, V1mag, b_eq));
    axis tight; yl = ylim; ylim([min(QvalsMVAr), max(QvalsMVAr)]);
    cb = colorbar; cb.Ticks=[0.25 0.75]; cb.TickLabels={'Feasible','Infeasible'};

    % ---- 6) Return parameters for traceability
    params = struct('Yred',Yred,'y_series',y_series,'b_eq',b_eq,'V1mag',V1mag);
end
