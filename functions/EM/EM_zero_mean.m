function [Ea, Eb, pa, pb, OPS, n_it] = EM_zero_mean(y, stop_criterion, n_it_max, Ea, Eb, Npa, Npb)

N = length(y);
abs_y2 = abs(y).^2; % 3*D

OPS = 0; % number of operations (detailed as comments next to each operation, e.g., 3*D+1)
n_it = 0;
stop = (n_it>=n_it_max);
while (~stop)
    n_it = n_it + 1;
    % expectation-maximization
    lika = exp(-abs_y2/Ea)*(Npa/Ea); % 3*D+1
    likb = exp(-abs_y2/Eb)*(Npb/Eb); % 3*D+1
    a = lika./(lika+likb); % 2*D
    b = 1-a; % 1*D      
    Npa = sum(a); % 1*D
    Npb = N-Npa; % 1
    Ea_new = (a'*abs_y2)/Npa; % 3*D
    Eb_new = (b'*abs_y2)/Npb; % 3*D
    OPS = OPS + 16*N + 3;
    
    % convergence
    param_change = abs((Ea-Ea_new)/Ea) + abs((Eb-Eb_new)/Eb);    
    stop = (param_change < stop_criterion)||(n_it>=n_it_max);
    OPS = OPS + 9;
    
    % update parameters
    Ea = Ea_new;
    Eb = Eb_new;
end
pa = Npa/N;
pb = Npb/N;
end
