%% simulator for blind parametric activity rate estimators (Fig. 3)
% Gaussian noise plus Gaussian signal with probability p
% 2022 (c) ag753@cornell.edu
% Edited on March 2023 by ag753@cornell.edu: clean up and comment the code
% for GitHub upload
% -----------------------------------------------------

clear all
addpath(genpath("functions")) % add functions to path

% rng(0) % set random seed for reproducibility
trials = 10000; % number of Monte-Carlo trials
D = 64; % dimension of the vectors (larger dimension gives less variance)
p_list = 0:0.01:0.5; % activity rate (low p means more sparse)
N0 = 1; % noise power (fixed for this simulation)
Es_list = [0.5, 10]; % signal power
q_list = [1 1 2 2]; % parameter for activity rate estimator
r_list = [2 Inf 4 Inf]; % parameter for activity rate estimator

% sample estimate (activity rate of each realization of the signal):
p_sample_list = NaN(trials,length(Es_list),length(p_list)); 
% blind parametric estimate:
p_qr_est_list = NaN(trials,length(Es_list),length(p_list),length(q_list)); 

for tt = 1:trials % Monte-Carlo trials
    for ee = 1:length(Es_list) % sweep signal power
        for pp = 1:length(p_list) % sweep activity rate
            p = p_list(pp);
                   
            n = sqrt(N0*0.5)*(randn(D,1)+1i*randn(D,1));
            s = sqrt(Es_list(ee)/p*0.5)*(randn(D,1)+1i*randn(D,1));
            sparse_idices = rand(D,1)>=p;
            s(sparse_idices) = 0;
            
            y = s + n;                                                     

            p_sample_list(tt,ee,pp) = 1-sum(sparse_idices)/D;    
            for qq = 1:length(q_list)
                p_qr_est_list(tt,ee,pp,qq) = estp(q_list(qq),r_list(qq),y,D);                 
            end                                                
        end
    end
end

%% Plots
show_variance = 0;
for SNR_idx = 1:2
    figure()
    if (show_variance)
        plot(p_list,p_list*nan);
        ColOrd = get(gca,'ColorOrder');
        plot(p_list,p_list,'k:','DisplayName','reference','LineWidth',2);
        hold on;
        shadedErrorBar(p_list,squeeze(p_sample_list(:,SNR_idx,:)),{@mean,@std},...
            'lineprops',{'color',ColOrd(1,:),'displayname','genie'},'transparent',true,'patchSaturation',0.1);
        for qq=1:length(q_list)
            shadedErrorBar(p_list,squeeze(p_qr_est_list(:,SNR_idx,:,qq)),{@mean,@std},...
                'lineprops',{'color',ColOrd(qq+1,:),...
                'displayname',['$\hat{p}(' num2str(q_list(qq)) ',' strrep(num2str(r_list(qq)),'Inf','\infty') ')$']},...
                'transparent',true,'patchSaturation',0.1);
        end
        legend('Interpreter','latex','EdgeColor','k')
        xlabel('$p=E\left[\|\mathbf{s}\|_0\right]/D$','Interpreter','latex');
        ylabel('estimate for $p$','Interpreter','latex');
        title(['SNR = ' num2str(Es_list(SNR_idx)/N0)]);
        grid on;
        hold off;
        ylim([0,0.5])
        set(gca,'FontSize',14);
    else
        plot(p_list,p_list,'k:','DisplayName','reference','LineWidth',2);
        hold on;
        for qq=1:length(q_list)
            plot(p_list,squeeze(mean(p_qr_est_list(:,SNR_idx,:,qq),1)),...
                'DisplayName',['$\hat{p}(' num2str(q_list(qq)) ',' strrep(num2str(r_list(qq)),'Inf','\infty') ')$'],...
                'LineWidth',2);
        end
        legend('Interpreter','latex','EdgeColor','k');
        xlabel('$p=E\left[\|\mathbf{s}\|_0\right]/D$',...
            'Color','k','Interpreter','latex');
        ylabel('estimate for $p$','Color','k', 'Interpreter','latex');
        title(['SNR = ' num2str(Es_list(SNR_idx)/N0)]);
        grid on;
        hold off;
        ylim([0,0.5])
        set(gca,'FontSize',14);
    end
end