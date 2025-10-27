% Q6 helper function
function feas = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, P7_vals, Q7_vals, toler, maxiter)
    % pf_feasibility_map - compute feasibility map for varying load at bus 7
    % Input:
    %   Y - admittance matrix
    %   is - slack bus index
    %   ipq, ipv - PQ and PV bus indices
    %   Pg, Qg - generation vectors
    %   Pd, Qd - demand vectors
    %   V0 - initial voltage magnitudes
    %   Sbase - base power (MVA)
    %   P7_vals, Q7_vals - load variation ranges for bus 7
    %   toler, maxiter - convergence parameters
    % Output:
    %   feas - feasibility matrix (boolean)

feas = false(numel(P7_vals), numel(Q7_vals));
for a = 1:numel(P7_vals)
    for b = 1:numel(Q7_vals)
        Pd_try = Pd; Qd_try = Qd;
        % Qd_try
        Pd_try(7) = P7_vals(a);
        Qd_try(7) = Q7_vals(b);
        % toler = 0.00001; maxiter = 10;
        % try
        [V, ~, ~, ~, N, ~] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd_try, Qd_try, V0, Sbase, toler, maxiter);
        if N < maxiter
            if all(abs(V)>=0.95 & abs(V)<=1.05)
                feas(a,b) = true;
            end
        end
        % catch
        %     feas(a,b) = false;
        % end
    end
end
end
% 
% function F = pf_feasibility_map(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, ...
%                                  Pd7_range, Qd7_range, toler, maxiter)
% 
% %--------------------------------------------------------------------------
% % Power Flow Feasibility Map at Bus 7
% %--------------------------------------------------------------------------
% % Inputs:
% %   Y          : System admittance matrix
% %   is, ipq, ipv : Bus type indicators
% %   Pg, Qg     : Generator active/reactive powers
% %   Pd, Qd     : Load active/reactive powers
% %   V0         : Initial voltage guess
% %   Sbase      : Base power (MVA)
% %   Pd7_range  : Range of active power at Bus 7 (MW)
% %   Qd7_range  : Range of reactive power at Bus 7 (Mvar)
% %   toler      : NR tolerance
% %   maxiter    : NR max iteration count
% %
% % Output:
% %   Generates a filled contour plot showing feasible region at Bus 7.
% %--------------------------------------------------------------------------
% 
% F = zeros(length(Pd7_range), length(Qd7_range));  % Feasibility matrix
% % toler = 0.00001; maxiter = 20;
% for a = 1:length(Pd7_range)
%     for k = 1:length(Qd7_range)
% 
%         Pd(7) = Pd7_range(a);
%         Qd(7) = Qd7_range(k);
% 
%         try
%             % Run power flow
%             [V_q6, delta_q6, Psl_q6, Qgv_q6, N_q6, time_q6] = ...
%                 nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);
% 
%             % Check convergence and voltage magnitude range
%             if N_q6 < maxiter
%                 Vmag = abs(V_q6);
%                 if all(Vmag >= 0.95 & Vmag <= 1.05)
%                     F(a,k) = 1;
%                 end
%             end
% 
%         catch
%             % If NR fails, mark as infeasible
%             F(a,k) = 0;
%         end
% 
%     end
% end
% 
