% Q3
function [V, delta, Psl, Qgv, N, time] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    tic; % Start timing the computation

    % Extract the dimensions of the system
    n = size(Y, 1);

    % Initialize voltage magnitudes and angles
    V = V0; % Initial voltage magnitudes
    delta = zeros(n, 1); % Initial voltage angles in radians

    % Calculate net power injections
    P = Pg - Pd; % Net active power injection (generation - demand)
    Q = Qg - Qd; % Net reactive power injection (generation - demand)

    % Initialize for iteration
    Vm = V; % Magnitude of voltage vector
    Theta = delta; % Angle of voltage vector in radians
    mismatchP = inf; % Initialize mismatch for active power to a large number
    mismatchQ = inf; % Initialize mismatch for reactive power to a large number
    N = 0; % Iteration counter

    % Precompute parts of the Jacobian that are constant
    Bp = -imag(Y(ipq, ipq));
    Bpp = -imag(Y(ipq, ipq));

    % Fast Decoupled iterations
    while (max(abs(mismatchP)) > toler || max(abs(mismatchQ)) > toler) && (N < maxiter)
        [mismatchP] = calcMismatchP(Y, Vm, Theta, P, ipq);
        dTheta = -Bp \ mismatchP;
        Theta(ipq) = Theta(ipq) + dTheta;
        
        [mismatchQ] = calcMismatchQ(Y, Vm, Theta, Q, ipq);
        dV = -Bpp \ mismatchQ;
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

    function [mismatchP] = calcMismatchP(Y, Vm, Theta, P, ipq)
        % Compute the active power mismatch
        I = Y * (Vm .* exp(1j * Theta));
        S = Vm .* conj(I);
        mismatchP = real(S(ipq)) - P(ipq);
    end

    function [mismatchQ] = calcMismatchQ(Y, Vm, Theta, Q, ipq)
        % Compute the reactive power mismatch
        I = Y * (Vm .* exp(1j * Theta));
        S = Vm .* conj(I);
        mismatchQ = imag(S(ipq)) - Q(ipq);
    end
end
