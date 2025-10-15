% Q1
function [V, delta, Psl, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
    tic; % start timing the computation
    
    % number of buses
    n = size(Y, 1);


    jp = [2:n]';
    psize = n - 1;
    jq = ipq;
    qsize = size(ipq, 1)
    jsize = psize + qsize; % the size N of the jacobian matrix
    
    % Initialize the complex voltages to the flat voltage profile
    profile = [zeros(1, psize), ones(1, qsize)]' % [delta | V]

    delta = abs(profile(1:psize));
    V = abs(profile(psize+1:end));

    % H = NaN(theta_size, theta_size);
    % % N = zeros(theta_size, pq_size);
    % % M = zeros(pq_size, theta_size);
    % L = NaN(pq_size, pq_size);
    % 
    % 
    % N = NaN(theta_size, theta_size);
    % M = NaN(theta_size, theta_size);

    % jac_size = size(approx, 1); 

    J = zeros(jsize, jsize);

    G = real(Y);
    B = imag(Y);
    
    % Net active and reactive power injection
    P = Pg - Pd; % Net active power injection (generation - demand)
    Q = Qg - Qd; % Net reactive power injection (generation - demand)

    
    Inj = [P(jp);Q(jq)];


    iter_count = 0; % Iteration counter



    while (iter_count < maxiter)


        delta_full = [0; delta];
    
        Vm = [V0; V]
        V_full = Vm .* exp(1i*delta_full);
    
        S = V_full .* conj(Y*V_full);
    
        
        mismatch = [P(jp)-real(S(jp));Q(jq)-imag(S(jq))]

        if (max(abs(mismatch)) < toler)
            break
        end

        % V = [V0; approx(psize+1:end)];
        % delta = [1; approx(1:psize)];
        dSdd = (-1i * diag(V_full) * conj(Y) * diag(conj(V_full)) + 1i * diag(conj(Y) * conj(V_full)) * diag(V_full));
        dSdVm = (diag(V_full) * conj(Y) * diag(exp(-1i * delta_full)) + diag(conj(Y) * conj(V_full)) * diag(exp(1i * delta_full)));
    
        H = real(dSdd);
        M = imag(dSdd);
        N_ = diag(Vm) * real(dSdVm); %rename N so it does not cause comflict with N
        L = diag(Vm) * imag(dSdVm);
    
        J(1:psize,1:psize) = H(jp,jp);
        J(1:psize,psize+1:end) = N_(jp,jq);
        J(psize+1:end,1:psize) = M(jq,jp);
        J(psize+1:end,psize+1:end) = L(jq,jq);

        if (max(abs(mismatch)) < toler)
            break
        end

        solution = linsolve(J,mismatch) %3 solve to correct


        delta = delta + solution(1:psize);
        V = V + V.*solution(psize+1:end);
        % profile = [profile(1:psize)+solution(1:psize); profile(psize+1:end)+profile(psize+1:end).*solution(psize+1:end)] %4 update profile
        
        % 
        % mismatch = Inj - approx; %new mismatch

        iter_count = iter_count+1

        J


    end


    % delta = abs(profile(1:psize));
    % V = abs(profile(psize+1:end));

    N = iter_count;

    Psl = delta; % not computed yet
    Qgv = V;

    time = toc;

    


  
end
