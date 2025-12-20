function results = dc_scopf(ifrom, ito, x, fmax, d, co, a, b, gmin, gmax, ngen, is)
% Security-constrained DC optimal power flow solver
% Uses quadratic programming for optimal dispatch with line constraints

    % Extract system dimensions
    num_buses = max([ifrom; ito]);
    num_lines = length(ifrom);
    num_gens = length(co);
    gen_locations = ngen(:);
    slack_bus = is;
    
    % Compute line admittances (susceptances)
    line_admittances = x .^ (-1);
    
    % Construct bus admittance matrix using sparse indexing
    Y_bus = construct_admittance_matrix(ifrom, ito, line_admittances, num_buses, num_lines);
    
    % Identify non-reference buses for reduced formulation
    active_buses = setdiff(1:num_buses, slack_bus);
    num_active = length(active_buses);
    Y_reduced = Y_bus(active_buses, active_buses);
    
    % Build power transfer distribution factor (PTDF) matrix
    incidence_mat = build_incidence_matrix(ifrom, ito, num_lines, num_buses);
    PTDF_mat = diag(line_admittances) * incidence_mat(:, active_buses);
    
    % Setup generator-to-bus mapping matrix
    gen_map = sparse(gen_locations, 1:num_gens, ones(num_gens,1), num_buses, num_gens);
    gen_map_reduced = gen_map(active_buses, :);
    demand_reduced = d(active_buses);
    
    % Formulate QP objective: minimize 0.5*x'*Q*x + c'*x
    Q_matrix = construct_cost_matrix(b, num_gens, num_active);
    c_vector = [a(:); zeros(num_active, 1)];
    
    % Construct equality constraints for power balance
    [A_eq, b_eq] = build_equality_constraints(gen_map_reduced, Y_reduced, ...
                                               demand_reduced, num_gens, num_active, d);
    
    % Construct inequality constraints for line flow limits  
    [A_ineq, b_ineq] = build_inequality_constraints(PTDF_mat, fmax, num_lines, num_gens);
    
    % Set bounds on decision variables (generation and angles)
    lower_bounds = [gmin(:); -inf(num_active, 1)];
    upper_bounds = [gmax(:); inf(num_active, 1)];
    
    % Solve quadratic program
    qp_options = optimset('Display', 'off');
    [solution, ~, flag, ~, multipliers] = quadprog(Q_matrix, c_vector, ...
        A_ineq, b_ineq, A_eq, b_eq, lower_bounds, upper_bounds, [], qp_options);
    
    if flag <= 0
        error('Optimization failed with exit flag: %d', flag);
    end
    
    % Parse solution vector
    generation = solution(1:num_gens);
    angles_reduced = solution(num_gens+1:end);
    
    % % Reconstruct full angle vector
    % angles_full = zeros(num_buses, 1);
    % angles_full(active_buses) = angles_reduced;
    % 
    % % Calculate line flows using PTDF
    line_flows = PTDF_mat * angles_reduced;
    
    % Compute total cost
    total_cost = sum(co(:) + a(:).*generation + 0.5*b(:).*(generation.^2));
    
    % Extract locational marginal prices from dual variables
    dual_equality = multipliers.eqlin;
    dual_buses = dual_equality(1:num_active);
    dual_balance = dual_equality(num_active + 1);
    
    lmp = zeros(num_buses, 1);
    lmp(active_buses) = dual_buses - dual_balance;
    lmp(slack_bus) = -dual_balance;
    
    % Calculate net power injections and merchandizing surplus
    net_injection = gen_map * generation - d;
    surplus = -sum(lmp .* net_injection);
    
    % Package results into output structure
    results.g = generation;
    % results.delta = angles_full;
    results.flow = line_flows;
    results.C = total_cost;
    results.LMP = lmp;
    results.MS = surplus;
    % results.lambda = multipliers;
    % results.pinj = net_injection;
    % results.ifrom = ifrom;
    % results.ito = ito;
    % results.refbus = slack_bus;
end

function Y = construct_admittance_matrix(from_bus, to_bus, admittances, n_bus, n_line)
    Y = zeros(n_bus, n_bus);
    for idx = 1:n_line
        bus_i = from_bus(idx);
        bus_j = to_bus(idx);
        y_val = admittances(idx);
        Y(bus_i, bus_i) = Y(bus_i, bus_i) + y_val;
        Y(bus_j, bus_j) = Y(bus_j, bus_j) + y_val;
        Y(bus_i, bus_j) = Y(bus_i, bus_j) - y_val;
        Y(bus_j, bus_i) = Y(bus_j, bus_i) - y_val;
    end
end

function A = build_incidence_matrix(from_bus, to_bus, n_line, n_bus)
    A = zeros(n_line, n_bus);
    for idx = 1:n_line
        A(idx, from_bus(idx)) = 1;
        A(idx, to_bus(idx)) = -1;
    end
end

function Q = construct_cost_matrix(cost_coeff, n_gen, n_angle)
    Q = blkdiag(diag(cost_coeff(:)), zeros(n_angle, n_angle));
end

function [A_eq, b_eq] = build_equality_constraints(G_red, Y_red, d_red, n_gen, n_active, d_full)
    A_nodal = [-G_red, Y_red];
    b_nodal = -d_red;
    A_global = [ones(1, n_gen), zeros(1, n_active)];
    b_global = sum(d_full);
    A_eq = [A_nodal; A_global];
    b_eq = [b_nodal; b_global];
end

function [A_ineq, b_ineq] = build_inequality_constraints(PTDF, f_lim, n_line, n_gen)
    A_ineq = [zeros(n_line, n_gen), PTDF; 
              zeros(n_line, n_gen), -PTDF];
    b_ineq = [f_lim(:); f_lim(:)];
end
