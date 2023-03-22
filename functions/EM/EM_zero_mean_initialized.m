function [Ea, Eb, pa, pb, OPS, n_it] = EM_zero_mean_initialized(y, stop_criterion, n_it_max)

D = length(y);

% initialize
Dpb = sum(abs(y))/max(abs(y)); % D * norm-1 / norm-infinity
Dpa = D-Dpb;
Ea = (y'*y)/D*0.4;
Eb = Ea+max((y'*y)/D-Ea,0)/(Dpb/D);

abs_y2 = abs(y).^2; % 3*D

OPS = 0;
n_it = 0;
stop = (n_it>=n_it_max);
while (~stop)
    n_it = n_it + 1;
    % expectation-maximization
    lika = exp(-abs_y2/Ea)*(Dpa/Ea); % 3*D+1
    likb = exp(-abs_y2/Eb)*(Dpb/Eb); % 3*D+1
    a = lika./(lika+likb); % 2*D
    b = 1-a; % 1*D      
    Dpa = sum(a); % 1*D
    Dpb = D-Dpa; % 1
    Ea_new = (a'*abs_y2)/Dpa; % 3*D
    Eb_new = (b'*abs_y2)/Dpb; % 3*D
    OPS = OPS + 16*D + 3;
    
    % convergence
    param_change = abs((Ea-Ea_new)/Ea) + abs((Eb-Eb_new)/Eb);
    stop = (param_change < stop_criterion)||(n_it>=n_it_max);
    OPS = OPS + 9;
    
    % update parameters
    Ea = Ea_new;
    Eb = Eb_new;
end
pa = Dpa/D;
pb = Dpb/D;
end
