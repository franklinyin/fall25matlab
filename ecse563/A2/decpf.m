% Q2
function [V, delta, Psl, Qgv, N, time] = decpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    % decpf - decoupled power flow (iterates with H and L blocks only)
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
    n = size(Y,1); G = real(Y); B = imag(Y);
    ipq = unique(ipq(:));
    ipv = setdiff(unique(ipv(:)), is);
    
    Pinj = (Pg - Pd)/Sbase;  Qinj = (Qg - Qd)/Sbase;
    
    V     = ones(n,1);
    delta = zeros(n,1);
    
    if numel(V0) == numel(ipv), V(ipv) = V0(:); end
    
    iang = [ipv; ipq]; iV = ipq;
    
    tic
    for N = 1:maxiter
        % compute P, Q
        P = zeros(n,1); Q = zeros(n,1);
        for i = 1:n
            for k = 1:n
                P(i) = P(i) + V(i)*V(k)*( G(i,k)*cos(delta(i)-delta(k)) + B(i,k)*sin(delta(i)-delta(k)) );
                Q(i) = Q(i) + V(i)*V(k)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
            end
        end
        dP = Pinj - P; dQ = Qinj - Q;
    
        if max( [norm(dP(iang),inf) norm(dQ(iV),inf)] ) < toler
            break
        end
    
        % build H and L only (as opposed to the full build in NR)
        np = numel(iang); nq = numel(iV);
        H = zeros(np,np); L = zeros(nq,nq);
        for a = 1:np
            i = iang(a);
            for b = 1:np
                k = iang(b);
                if i==k, H(a,b) = -Q(i) - V(i)^2*B(i,i);
                else,    H(a,b) = V(i)*V(k)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
                end
            end
        end
        for a = 1:nq
            i = iV(a);
            for b = 1:nq
                k = iV(b);
                if i==k, L(a,b) =  Q(i)/V(i) - B(i,i)*V(i);
                else,    L(a,b) =  V(i)*( G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)) );
                end
            end
        end
    
        % solve decoupled systems
        ddelta = H \ dP(iang);
        dV     = L \ dQ(iV);
    
        delta(iang) = delta(iang) + ddelta;
        V(iV)       = V(iV)       + dV;
        delta(is)   = 0;
        if numel(V0) == numel(ipv), V(ipv) = V0(:); end
    end
    time = toc;
    
    % Psl and Qgv as in NR
    Psl = ( P(is) + (Pd(is)-Pg(is))/Sbase ) * Sbase;
    Qgv = zeros(n,1);
    for k = 1:numel(ipv)
        i = ipv(k); Qgv(i) = ( Q(i) + Qd(i)/Sbase ) * Sbase;
    end
end
