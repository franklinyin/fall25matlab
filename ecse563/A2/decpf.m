% Q2
function [V, delta, Psl, Qgv, N, time] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    tic; % Start timing the computation
    
    % Number of buses
    n = length(Y);
    
    % Initialize
    V = V0; % Voltage magnitudes
    delta = zeros(n, 1); % Voltage angles in radians
    P = Pg - Pd; % Net active power injection (generation - demand)
    Q = Qg - Qd; % Net reactive power injection (generation - demand)
    
    % Initialize for iteration
    Vm = V; % Magnitude of voltage vector
    Theta = delta; % Angle of voltage vector in radians
    mismatchP = inf; % Initialize mismatch for active power to a large number
    mismatchQ = inf; % Initialize mismatch for reactive power to a large number
    N = 0; % Iteration counter
    
    % Decoupled iteration
    while (max(abs(mismatchP)) > toler || max(abs(mismatchQ)) > toler) && (N < maxiter)
        [mismatchP, Jp] = calcMismatchP(Y, Vm, Theta, P, ipq, is);
        dTheta = -Jp \ mismatchP;
        Theta(ipq) = Theta(ipq) + dTheta;
        
        [mismatchQ, Jq] = calcMismatchQ(Y, Vm, Theta, Q, ipq);
        dV = -Jq \ mismatchQ;
        Vm(ipq) = Vm(ipq) + dV;
        
        % Increment iteration counter
        N = N + 1;
    end
    
    % Outputs after convergence or reaching max iterations
    V = Vm;
    delta = Theta;
    Psl = P(is) - sum(Y(is,:) .* (Vm .* exp(1j * Theta))); % Slack bus active power
    Qgv = Q(ipv); % Reactive power at PV nodes
    
    % Timing
    time = toc;

    function [mismatchP, Jp] = calcMismatchP(Y, Vm, Theta, P, ipq, is)
        % Compute the active power mismatch and the Jacobian matrix for Theta
        I = Y * (Vm .* exp(1j * Theta));
        S = Vm .* conj(I);
        mismatchP = real(S(ipq)) - P(ipq);
        Jp = -diag(Vm(ipq)) * imag(Y(ipq, ipq)) * diag(cos(Theta(ipq))) - ...
             diag(Vm(ipq)) * real(Y(ipq, ipq)) * diag(sin(Theta(ipq)));
    end

    function [mismatchQ, Jq] = calcMismatchQ(Y, Vm, Theta, Q, ipq)
        % Compute the reactive power mismatch and the Jacobian matrix for Vm
        I = Y * (Vm .* exp(1j * Theta));
        S = Vm .* conj(I);
        mismatchQ = imag(S(ipq)) - Q(ipq);
        Jq = diag(real(Y(ipq, ipq))) * diag(Vm(ipq));
    end
end
