% Q3
function [V, delta, Psl, Qgv, N, time] = fastdecpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    % fastdecpf - Stott's fast-decoupled power flow (B' and B'' constant)
    % Input:
    %   Y - admittance matrix
    %   is - slack bus index
    %   ipq, ipv - PQ and PV bus indices
    %   Pg, Qg - generation vectors
    %   Pd, Qd - demand vectors
    %   V0 - initial voltage magnitudes
    %   Sbase - base power (MVA)
    %   toler, maxiter - convergence parameters
    % Output:
    %   V - voltage magnitudes (p.u.)
    %   delta - voltage angles (radians)
    %   Psl - slack bus active power (MW)
    %   Qgv - PV bus reactive power (MVAr)
    %   N - number of iterations
    %   time - computation time (seconds)
   
    %
    n = size(Y,1); B = imag(Y);
    ipq = unique(ipq(:));
    ipv = setdiff(unique(ipv(:)), is);
    
    Pinj = (Pg - Pd)/Sbase;  Qinj = (Qg - Qd)/Sbase;
    
    % build B' from series susceptances only: diag = sum(|offdiag|), offdiag = +b_series
    Bseries = zeros(n);
    for i = 1:n
        for k = 1:n
            if i~=k, Bseries(i,k) = B(i,k); end % off-diag of Y is series-only
        end
        Bseries(i,i) = -sum(Bseries(i,[1:i-1 i+1:n])); % make Laplacian
    end
    Bprime = -Bseries; % conventional B'
    Bpp    = -B; % conventional B'' uses full -imag(Y)
    
    % State
    V     = ones(n,1);
    delta = zeros(n,1);
    if numel(V0)==numel(ipv), V(ipv) = V0(:); end
    
    iang = [ipv; ipq];  % non-slack angles
    iV   = ipq;
    
    % Reduced matrices
    Bp  = Bprime(iang, iang);
    Bpp = Bpp(iV, iV);
    
    tic
    for N = 1:maxiter
        % compute P, Q
        P = zeros(n,1); Q = zeros(n,1);
        for i = 1:n
            for k = 1:n
                Gik = real(Y(i,k));
                Bik = imag(Y(i,k));
                P(i) = P(i) + V(i)*V(k)*( Gik*cos(delta(i)-delta(k)) + Bik*sin(delta(i)-delta(k)) );
                Q(i) = Q(i) + V(i)*V(k)*( Gik*sin(delta(i)-delta(k)) - Bik*cos(delta(i)-delta(k)) );
            end
        end
        dP = Pinj - P; dQ = Qinj - Q;
    
        if max( [norm(dP(iang),inf) norm(dQ(iV),inf)] ) < toler
            break
        end
    
        % Normalize RHS by |V|
        rhsP = dP(iang) ./ V(iang);
        rhsQ = dQ(iV)   ./ V(iV);
    
        % Solve with constant B' and B''
        ddelta = Bp  \ rhsP;
        dV     = Bpp \ rhsQ;
    
        delta(iang) = delta(iang) + ddelta;
        V(iV)       = V(iV)       + dV;
    
        delta(is) = 0;
        if numel(V0)==numel(ipv), V(ipv) = V0(:); end
    end
    time = toc;
    
    % Psl and Qgv as in NR
    Psl = ( P(is) + (Pd(is)-Pg(is))/Sbase ) * Sbase;
    Qgv = zeros(n,1);
    for k = 1:numel(ipv)
        i = ipv(k); Qgv(i) = ( Q(i) + Qd(i)/Sbase ) * Sbase;
    end
end
