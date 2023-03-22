function sim_1_2_4_and_5(trials)
%% Simulator for blind estimators: Figures 1, 2, 4 and 5.
% 2020 (c) studer@ethz.ch and ag753@cornell.edu
% Edited on March 2023 by ag753@cornell.edu: clean up and comment the code
% for GitHub upload
% -----------------------------------------------------
addpath(genpath("functions")) % add functions to path

par.est_type_list = {'genie','blind_nonparametric','EM','blind_parametric_LB_avg_UB','accelerated_EM','baseline_EM'};

for bb = 1:length(par.est_type_list)
    par.est_type = par.est_type_list{bb};

    rng(0) % set random seed for reproducibility
    par.trials = trials; % number of Monte-Carlo trials
    par.tt_realization = 1; % which realization to save for the MSE example
    par.D_list = [64 256]; % dimension (larger dimension gives less variance)
    par.p_list = [0.1,0.4]; % list of activity rate values (fraction of nonzeros)
    par.N0 = 1; % noise power (fixed for this simulation)
    par.Es_list = linspace(0,10,101); % list of signal power values

    % predefined settings for each estimator type
    switch par.est_type
        case {'genie','blind_nonparametric'}
            par.q_list = 1;
            par.r_list = 1;
            par.EM.stop_criterion = 0;
            par.EM.n_it_max_list = 1;
        case {'blind_parametric_LB_avg_UB'}
            par.q_list = [NaN 1]; % par.q_list = [NaN NaN 1 2];
            par.r_list = [0.25 Inf]; % par.r_list = [0 0.25 Inf 4];
            par.EM.stop_criterion = 0;
            par.EM.n_it_max_list = 1;
        case {'EM'}
            par.q_list = 1;
            par.r_list = 1;
            par.EM.stop_criterion = 0.001;
            par.EM.n_it_max_list = 30;
        case {'accelerated_EM','baseline_EM'}
            par.q_list = [NaN 1]; % par.q_list = [NaN NaN 1 2];
            par.r_list = [0.25 Inf]; % par.r_list = [0 0.25 Inf 4];
            par.EM.stop_criterion = 0;
            par.EM.n_it_max_list = 2.^(0:7);
            par.Es_list = par.Es_list([2,51,101]); % only simulate 0.1 and 5
            par.p_list = [0.1,0.4]; % sparsity
            par.D_list = 256;
    end

    par.SNR_list = par.Es_list/par.N0;

    complexity = zeros(par.trials,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list)); % number of operations
    n_it = zeros(par.trials,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list)); % number of iterations of the EM algorithm

    % initialize
    N0_est_list = NaN(par.trials,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    Es_est_list = NaN(par.trials,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    SNR_est_list = NaN(par.trials,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));

    % one instance of D-dimensional vectors tau and MSE per list of parameters
    tau_est_list = NaN(1,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list),max(par.D_list));
    MSE_est_list = NaN(1,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list),max(par.D_list));

    % GO
    for tt = 1:par.trials
        for dd = 1:length(par.D_list)
            D = par.D_list(dd);
            for pp = 1:length(par.p_list)
                p = par.p_list(pp);
                for ee = 1:length(par.Es_list)
                    Es = par.Es_list(ee);

                    % generate Gaussian noise
                    n = sqrt(par.N0*0.5)*(randn(D,1)+1i*randn(D,1));

                    % generate sparse signal
                    s = sqrt(Es/p*0.5)*(randn(D,1)+1i*randn(D,1));
                    sparse_idices = rand(D,1)>=p;
                    s(sparse_idices) = 0;

                    %noisy signal
                    y = s+n;

                    rand; rand; % not important, only needed to reproduce the exact realizations of fig.  in the paper.

                    for qq=1:length(par.q_list)
                        q = par.q_list(qq);
                        r = par.r_list(qq);

                        % estimate the sparsity
                        if(isnan(q)) % fix p_est = r
                            if (r>0 && r<1)
                                p_est = r;
                            else
                                error('sparsity estimator not defined');
                            end
                        else % p(q,r) with vector norms
                            p_est = min(0.499,estp(q,r,y,D));
                        end


                        for ii = 1:length(par.EM.n_it_max_list) % iterate to study the convergence of EM
                            n_it_max = par.EM.n_it_max_list(ii);

                            % estimate the noise and signal power
                            switch par.est_type
                                case 'blind_nonparametric' % nonparametric
                                    [N0estfloor,N0estceil,complexity(tt,ee,pp,dd,qq,ii)] = ...
                                        qselect_modified(abs(y).^2/log(2),floor((D+1)/2),ceil((D+1)/2));
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = (N0estfloor + N0estceil)/2; % same as median(abs(y).^2/log(2))
                                case 'EM'
                                    [Ea,Eb,pa,pb,complexity(tt,ee,pp,dd,qq,ii),n_it(tt,ee,pp,dd,qq,ii)] = ...
                                        EM_zero_mean_initialized(y,par.EM.stop_criterion,n_it_max);
                                    EMvariances = [Ea Eb];
                                    EMprobabilities = [pa pb];
                                    [N0est, N0_idx] = min(EMvariances);
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = N0est;
                                    Es_idx = 3 - N0_idx;
                                    Esest = (EMvariances(Es_idx)-EMvariances(N0_idx))*...
                                        EMprobabilities(Es_idx);
                                case 'accelerated_EM'
                                    N0_median = median(abs(y).^2)/log(2);
                                    Es_median = max(0,(y'*y)/(D)-N0_median);
                                    % initialization for EM
                                    Ea = N0_median;
                                    Npb = p_est*D;
                                    Eb = Ea+Es_median/p_est;
                                    Npa = D-Npb;
                                    [Ea,Eb,pa,pb,complexity(tt,ee,pp,dd,qq,ii),n_it(tt,ee,pp,dd,qq,ii)] = ...
                                        EM_zero_mean(y,par.EM.stop_criterion,n_it_max,Ea,Eb,Npa,Npb);
                                    EMvariances = [Ea Eb];
                                    EMprobabilities = [pa pb];
                                    [N0est, N0_idx] = min(EMvariances);
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = N0est;
                                    Es_idx = 3 - N0_idx;
                                    Esest = (EMvariances(Es_idx)-EMvariances(N0_idx))*...
                                        EMprobabilities(Es_idx);
                                case 'baseline_EM'
                                    N0fix = (y'*y)/D/6; % N0 by assumung SNR=5 (or ~7dB)
                                    Esfix = max(0,(y'*y)/D-N0fix);
                                    % initialization for EM
                                    Ea = N0fix;
                                    Npb = p_est*D;
                                    Eb = Ea+Esfix/p_est;
                                    Npa = D-Npb;
                                    [Ea,Eb,pa,pb,complexity(tt,ee,pp,dd,qq,ii),n_it(tt,ee,pp,dd,qq,ii)] = ...
                                        EM_zero_mean(y,par.EM.stop_criterion,n_it_max,Ea,Eb,Npa,Npb);
                                    EMvariances = [Ea Eb];
                                    EMprobabilities = [pa pb];
                                    [N0est, N0_idx] = min(EMvariances);
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = N0est;
                                    Es_idx = 3 - N0_idx;
                                    Esest = (EMvariances(Es_idx)-EMvariances(N0_idx))*...
                                        EMprobabilities(Es_idx);
                                case 'blind_parametric_LB_avg_UB'
                                    N0_median = median(abs(y).^2)/log(2);
                                    SNR_median = max(0,(y'*y)/(D*N0_median)-1);
                                    N0LB = N0_median/min(log((2-2*p_est)/(1-2*p_est))/log(2),(1+SNR_median));
                                    N0UB = N0_median*(1-p_est+p_est^2/(p_est+SNR_median));
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = (N0LB+N0UB)/2;
                                case 'genie'
                                    N0_est_list(tt,ee,pp,dd,qq,ii) = mean(abs(n).^2);
                                    Esest = mean(abs(s).^2);
                                otherwise
                                    error('est_type not defined')
                            end

                            if exist('Esest','var')
                                Es_est_list(tt,ee,pp,dd,qq,ii) = max(Esest,0);
                            else
                                Es_est_list(tt,ee,pp,dd,qq,ii) = max(norm(y,2)^2/D-N0_est_list(tt,ee,pp,dd,qq,ii),0);
                            end
                            clear Esest
                            SNR_est_list(tt,ee,pp,dd,qq,ii) = Es_est_list(tt,ee,pp,dd,qq,ii)./N0_est_list(tt,ee,pp,dd,qq,ii);

                            % compute MSE
                            if(tt==par.tt_realization)
                                switch par.est_type
                                    case 'genie'
                                        tau_est_list(1,ee,pp,dd,qq,ii,1:D) = sort(abs(y),'ascend');
                                        for aa=1:(D-1)
                                            eta_y_tau = (y(y~=0)./abs(y(y~=0))).*...
                                                max(abs(y(y~=0))-tau_est_list(1,ee,pp,dd,qq,ii,aa),0);
                                            MSE_est_list(1,ee,pp,dd,qq,ii,aa) = (s-eta_y_tau)'*(s-eta_y_tau)/D;
                                        end
                                    otherwise
                                        [~,~,~,tau_est_list(1,ee,pp,dd,qq,ii,1:D),MSE_est_list(1,ee,pp,dd,qq,ii,1:D)] = ...
                                            BEACHES(y,N0_est_list(tt,ee,pp,dd,qq,ii),'none');
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %% save only mean and standard deviation

    N0_est_list_mv = NaN(2,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    N0_est_list_mv(1,:,:,:,:,:) = mean(N0_est_list,1);
    N0_est_list_mv(2,:,:,:,:,:) = std(N0_est_list,0,1);

    N0_percent_error = abs(N0_est_list-par.N0)./par.N0;
    N0_percent_error_mv = NaN(2,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    N0_percent_error_mv(1,:,:,:,:,:) = mean(N0_percent_error,1);
    N0_percent_error_mv(2,:,:,:,:,:) = std(N0_percent_error,0,1);

    clear N0_percent_error

    clear N0_est_list

    Es_est_list_mv = NaN(2,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    Es_est_list_mv(1,:,:,:,:,:) = mean(Es_est_list,1);
    Es_est_list_mv(2,:,:,:,:,:) = std(Es_est_list,0,1);

    clear Es_est_list

    SNR_est_list_mv = NaN(2,length(par.Es_list),length(par.p_list),...
        length(par.D_list),length(par.q_list),length(par.EM.n_it_max_list));
    SNR_est_list_mv(1,:,:,:,:,:) = mean(SNR_est_list,1);
    SNR_est_list_mv(2,:,:,:,:,:) = std(SNR_est_list,0,1);

    clear SNR_est_list

    save(['./results/' par.est_type '_' num2str(par.trials) '.mat'],...
        'par','N0_est_list_mv','Es_est_list_mv','SNR_est_list_mv','tau_est_list',...
        'MSE_est_list','N0_percent_error_mv','-v7.3');

end
end