%% Script to generate Figures 1, 2, 4 and 5.
% 2020 (c) ag753@cornell.edu
% Edited on March 2023 by ag753@cornell.edu: clean up and comment the code
% for GitHub upload
% -----------------------------------------------------

clear all
%% select simulation options

est_list = {'blind_nonparametric','blind_parametric_LB_avg_UB','EM','genie'}; % list of estimators to be plotted
trials = 10000; % number of Monte-Carlo trials

show_variance = 1; % if show_variance = 0, only show mean performance

addpath(genpath('./functions'));

%% run simulation and save results
generate_results = 1; % can be set to zero if data is already generated and one only wants to plot
if generate_results
    sim_1_2_4_and_5(trials);
end

%% plotting options (colors, linestyles, labels, etc.)

pltopts = my_plotting_options_synthetic; % load a bunch of plotting options from separate file
my_interpreter = 'latex';

%% noise power
show_qq_indices = 2:2;
for pp=1:2
    for dd=1:2
        % noise power
        figure()
        for kk=1:length(est_list)
            if (show_variance ||~strcmp(est_list{kk},'genie'))
                my_variables = {'par','N0_est_list_mv'};
                load(['./results/' est_list{kk} '_' num2str(trials) '.mat'],my_variables{:});
                for qq=1:length(par.q_list)
                    if length(par.q_list)>1
                        if isnan(par.q_list(qq))
                            if (par.r_list(qq)>0 && par.r_list(qq)<1)
                                mylegend = [pltopts.est_label_N0{kk}...
                                    'p=' num2str(par.r_list(qq)) ')$'];
                            end
                        else
                            mylegend = [pltopts.est_label_N0{kk}...
                                num2str(par.q_list(qq)) ','...
                                strrep(num2str(par.r_list(qq)),'Inf','\infty') '))$'];
                        end
                    else
                        mylegend = pltopts.est_label_N0{kk};
                    end
                    if (length(par.q_list)==1 || sum(qq == show_qq_indices))
                        if (show_variance)
                            s = shadedErrorBar(par.SNR_list,squeeze(N0_est_list_mv(1,:,pp,dd,qq,1)),squeeze(N0_est_list_mv(2,:,pp,dd,qq,1)),...
                                'lineprops',{'color', pltopts.color(est_list{kk}),...
                                'displayname',mylegend,...
                                'linestyle',pltopts.linestyle(est_list{kk})},...
                                'transparent',true,'patchSaturation',0.1);
                            set(s.edge,'LineWidth',1,'LineStyle',pltopts.linestyle(est_list{kk}))
                            s.mainLine.LineWidth = 2;
                            s.mainLine.Marker = pltopts.marker(est_list{kk});
                            s.mainLine.MarkerIndices= 2:(length(par.Es_list)-1)/10:(length(par.Es_list));
                            hold on
                        else
                            plot(par.SNR_list,squeeze(N0_est_list_mv(1,:,pp,dd,qq,end)),...
                                'Color', pltopts.color(est_list{kk}),...
                                'LineStyle',pltopts.linestyle(est_list{kk}),... % pltopts.linestyle(est_list{kk}),...
                                'Marker',pltopts.marker(est_list{kk}),...
                                'DisplayName',mylegend,...
                                'MarkerIndices',2:(length(par.Es_list)-1)/10:(length(par.Es_list)),...
                                'LineWidth',2);
                            hold on
                        end
                    end
                end
                clear('N0_est_list');
            end
        end
        plot(par.SNR_list,par.N0*ones(length(par.Es_list),1),'k:',...
            'DisplayName',pltopts.est_label_N0{end},'LineWidth',2);
        hold off
        grid on
        ylim([0.5 3])

        ylabel('estimate for $N_0$','Color','k','Interpreter',my_interpreter);
        legend('Interpreter',my_interpreter,'Location','northwest','EdgeColor','k')
        box on;
        xlabel('$\textit{SNR} = E_s/N_0$','Color','k','Interpreter',my_interpreter);
        title(['D=' num2str(par.D_list(dd)) ', p=' num2str(par.p_list(pp))]);
        set(gca,'FontSize',14);
    end
end

%% signal power and SNR
show_qq_indices = 2:2;
for pp=1:1
    for dd=1:1
        % signal power
        figure()
        for kk=1:length(est_list)
            if (show_variance ||~strcmp(est_list{kk},'genie'))
                my_variables = {'par','Es_est_list_mv'};
                load(['./results/' est_list{kk} '_' num2str(trials) '.mat'],my_variables{:});
                for qq=1:length(par.q_list)
                    if length(par.q_list)>1
                        if isnan(par.q_list(qq))
                            if (par.r_list(qq)>0 && par.r_list(qq)<1)
                                mylegend = [pltopts.est_label_Es{kk}...
                                    'p=' num2str(par.r_list(qq)) ')$'];
                            else
                                p_est_type = {'\sim U(0,1/2)','=tanh','=BEACHES'};
                                mylegend = [pltopts.est_label_Es{kk}...
                                    'p' p_est_type{par.r_list(qq)+1} ')$'];
                            end
                        else
                            mylegend = [pltopts.est_label_Es{kk}...
                                num2str(par.q_list(qq)) ','...
                                strrep(num2str(par.r_list(qq)),'Inf','\infty') '))$'];
                        end
                    else
                        mylegend = pltopts.est_label_Es{kk};
                    end
                    if (length(par.q_list)==1 || sum(qq == show_qq_indices))
                        if (show_variance)
                            s = shadedErrorBar(par.SNR_list,squeeze(Es_est_list_mv(1,:,pp,dd,qq,1)),squeeze(Es_est_list_mv(2,:,pp,dd,qq,1)),...
                                'lineprops',{'color', pltopts.color(est_list{kk}),...
                                'displayname',mylegend,...
                                'linestyle',pltopts.linestyle(est_list{kk})},...
                                'transparent',true,'patchSaturation',0.1);
                            set(s.edge,'LineWidth',1,'LineStyle',pltopts.linestyle(est_list{kk}))
                            s.mainLine.LineWidth = 2;
                            s.mainLine.Marker = pltopts.marker(est_list{kk});
                            s.mainLine.MarkerIndices = 2:(length(par.Es_list)-1)/10:(length(par.Es_list));
                            hold on
                        else
                            plot(par.SNR_list,squeeze(Es_est_list_mv(1,:,pp,dd,qq,1)),...
                                'Color', pltopts.color(est_list{kk}),...
                                'LineStyle',pltopts.linestyle(est_list{kk}),...
                                'DisplayName',mylegend,...
                                'Marker',pltopts.marker(est_list{kk}),...
                                'MarkerIndices',2:(length(par.Es_list)-1)/10:(length(par.Es_list)),...
                                'LineWidth',2);
                            hold on
                        end
                    end
                end
                clear('Es_est_list');
            end
        end
        plot(par.SNR_list,par.Es_list,'k:',...
            'DisplayName',pltopts.est_label_Es{end},'LineWidth',2);
        hold off
        grid on
        ylim([0 max(par.Es_list)*(show_variance*0.5+1)])
        ylabel('estimate for $E_s$','Color','k','Interpreter',my_interpreter);
        legend('Interpreter',my_interpreter,'Location','northwest','EdgeColor','k')
        box on;
        xlabel('$\textit{SNR} = E_s/N_0$','Color','k','Interpreter',my_interpreter);
        title(['D=' num2str(par.D_list(dd)) ', p=' num2str(par.p_list(pp))]);
        set(gca,'FontSize',14);

        % SNR
        figure()
        for kk=1:length(est_list)
            if (show_variance ||~strcmp(est_list{kk},'genie'))
                my_variables = {'par','SNR_est_list_mv'};
                load(['./results/' est_list{kk} '_' num2str(trials) '.mat'],my_variables{:});
                for qq=1:length(par.q_list)
                    if length(par.q_list)>1
                        if isnan(par.q_list(qq))
                            if (par.r_list(qq)>0 && par.r_list(qq)<1)
                                mylegend = [pltopts.est_label_SNR{kk}...
                                    'p=' num2str(par.r_list(qq)) ')$'];
                            else
                                p_est_type = {'\sim U(0,1/2)','=tanh','=BEACHES'};
                                mylegend = [pltopts.est_label_SNR{kk}...
                                    'p' p_est_type{par.r_list(qq)+1} ')$'];
                            end
                        else
                            mylegend = [pltopts.est_label_SNR{kk}...
                                num2str(par.q_list(qq)) ','...
                                strrep(num2str(par.r_list(qq)),'Inf','\infty') '))$'];
                        end
                    else
                        mylegend = pltopts.est_label_SNR{kk};
                    end
                    if (length(par.q_list)==1 || sum(qq == show_qq_indices))
                        if (show_variance)
                            s = shadedErrorBar(par.SNR_list,squeeze(SNR_est_list_mv(1,:,pp,dd,qq,1)),squeeze(SNR_est_list_mv(2,:,pp,dd,qq,1)),...
                                'lineprops',{'color', pltopts.color(est_list{kk}),...
                                'linestyle', pltopts.linestyle(est_list{kk}),...
                                'displayname',mylegend},...
                                'transparent',true,'patchSaturation',0.1);
                            set(s.edge,'LineWidth',1,'LineStyle',pltopts.linestyle(est_list{kk}))
                            s.mainLine.LineWidth = 2;
                            s.mainLine.Marker = pltopts.marker(est_list{kk});
                            s.mainLine.MarkerIndices = 2:(length(par.Es_list)-1)/10:(length(par.Es_list));
                            hold on
                        else
                            plot(par.SNR_list,squeeze(SNR_est_list_mv(1,:,pp,dd,qq,1)),...
                                'Color', pltopts.color(est_list{kk}),...
                                'LineStyle',pltopts.linestyle(est_list{kk}),...
                                'DisplayName',mylegend,...
                                'Marker',pltopts.marker(est_list{kk}),...
                                'MarkerIndices',2:(length(par.Es_list)-1)/10:(length(par.Es_list)),...
                                'LineWidth',2);
                            hold on
                        end
                    end
                end
                clear('SNR_est_list');
            end
        end
        plot(par.SNR_list,par.SNR_list,'k:',...
            'DisplayName',pltopts.est_label_SNR{end},'LineWidth',2);
        hold off
        grid on
        ylim([0 max(par.SNR_list)*(show_variance*0.5+1)])
        ylabel('estimate for $\textit{SNR}$','Color','k','Interpreter',my_interpreter);
        legend('Interpreter',my_interpreter,'Location','northwest','EdgeColor','k')
        box on;
        xlabel('$\textit{SNR} = E_s/N_0$','Color','k','Interpreter',my_interpreter);
        title(['D=' num2str(par.D_list(dd)) ', p=' num2str(par.p_list(pp))]);
        set(gca,'FontSize',14);
    end
end

%% MSE vs tau for a given realization
show_qq_indices = 2:2;
tt = 1;
for ss=[6,101]
    for pp=2:2
        for dd=2:2
            figure()
            for kk=1:length(est_list)
                my_variables = {'par','MSE_est_list','tau_est_list'};
                load(['./results/' est_list{kk} '_' num2str(trials) '.mat'],my_variables{:});
                for qq=1:length(par.q_list)
                    if(length(par.q_list)>1)
                        mylegend = [pltopts.est_label_MSE{kk}...
                            num2str(par.q_list(qq)) ','...
                            strrep(num2str(par.r_list(qq)),'Inf','\infty') '))$'];
                    else
                        mylegend = pltopts.est_label_MSE{kk};
                    end
                    if (length(par.q_list)==1 || sum(qq == show_qq_indices))
                        tau = squeeze(tau_est_list(tt,ss,pp,dd,qq,1,:));
                        MSE = squeeze(MSE_est_list(tt,ss,pp,dd,qq,1,:));
                        [~,mymarkerindices,~] = unique(round(tau*15/round(max(tau)))/(15/round(max(tau))));
                        mymarkerindices(1) = mymarkerindices(1)+1;
                        mymarkerindices(2) = mymarkerindices(2)+4;
                        plot(tau,MSE,...
                            'LineWidth',2,...
                            'DisplayName',mylegend,...
                            'Color', pltopts.color(est_list{kk}),...
                            'LineStyle',pltopts.linestyle(est_list{kk}),...
                            'Marker',pltopts.marker(est_list{kk}),...
                            'MarkerIndices',mymarkerindices);
                        hold on;
                    end
                end
                clear('MSE_est_list','tau_est_list');
            end
            legend('Interpreter',my_interpreter,'EdgeColor','k','Location','northwest');
            xlabel('$\tau$','Color','k','Interpreter',my_interpreter);
            ylabel('estimate for $\textit{MSE}$','Color','k','Interpreter',my_interpreter);
            xmin = tau(1);
            xmax = tau(end-1);
            xlim([xmin,xmax]);
            grid on
            title(['SNR=' num2str(par.SNR_list(ss)) ', D=' num2str(par.D_list(dd)) ', p=' num2str(par.p_list(pp))]);
            set(gca,'FontSize',14);
        end
    end
end


%% accelerated EM
est_list = {'baseline_EM','accelerated_EM'};
for ss=[1,2] % [2,51,101]
    for pp=[1,2]
        for dd=1:1
            % error in the noise power
            figure()
            for kk=1:length(est_list)
                if (show_variance ||~strcmp(est_list{kk},'genie'))
                    my_variables = {'par','N0_percent_error_mv'};
                    load(['./results/' est_list{kk} '_' num2str(trials) '.mat'],my_variables{:});
                    % for qq=1:length(par.q_list)
                    for qq=[1,2]
                        if length(par.q_list)>1
                            if isnan(par.q_list(qq))
                                p_est_type = {'\sim U(0,1/2)','tanh','BEACHES'};
                                if (par.r_list(qq)>0 && par.r_list(qq)<1)
                                    mydisplayname = [pltopts.display_name(est_list{kk})...
                                        '=0.25$'];
                                else
                                    mydisplayname = [pltopts.display_name(est_list{kk})...
                                        p_est_type{par.r_list(qq)+1} '$'];
                                end
                            else
                                mydisplayname = [pltopts.display_name(est_list{kk})...
                                    '=\hat{p}(' num2str(par.q_list(qq)) ','...
                                    strrep(num2str(par.r_list(qq)),'Inf','\infty') ')$'];
                            end
                        else
                            mydisplayname = pltopts.est_label_N0{kk};
                        end
                        if (show_variance)
                            s = shadedErrorBar(log2(par.EM.n_it_max_list),squeeze(N0_percent_error_mv(1,ss,pp,dd,qq,:)),squeeze(N0_percent_error_mv(2,ss,pp,dd,qq,:)),...
                                'lineprops',{'color', pltopts.color(est_list{kk}),...
                                'displayname',mydisplayname,...
                                'linestyle',pltopts.mylinestyle2list{qq}},...
                                'transparent',true,'patchSaturation',0.1);
                            set(s.edge,'LineWidth',1,'LineStyle',pltopts.mylinestyle2list{qq})
                            s.mainLine.LineWidth = 2;
                            s.mainLine.Marker = pltopts.marker(est_list{kk});
                            hold on
                        else
                            plot(log2(par.EM.n_it_max_list),squeeze(N0_percent_error_mv(1,ss,pp,dd,qq,:)),...
                                'Color', pltopts.color(est_list{kk}),...
                                'LineStyle',pltopts.mylinestyle2list{qq},...
                                'Marker',pltopts.marker(est_list{kk}),...
                                'DisplayName',mydisplayname,...
                                'LineWidth',2);
                            hold on
                        end
                    end
                    clear('N0_est_list');
                end
            end
            hold off
            grid on
            ylim([0 0.8])
            xlabel('number of iterations $K$','Color','k','Interpreter',my_interpreter);
            ylabel('relative error $\varepsilon^\textrm{EM}$','Color','k','Interpreter',my_interpreter);
            legend('Interpreter',my_interpreter,'EdgeColor','k')
            box on;
            set (gca, 'XTick', log2(par.EM.n_it_max_list));
            set (gca, 'XTickLabel', par.EM.n_it_max_list);
            xlim([0 max(log2(par.EM.n_it_max_list))]);
            title(['SNR=' num2str(par.SNR_list(ss)) ', D=' num2str(par.D_list(dd)) ', p=' num2str(par.p_list(pp))]);
            set(gca,'FontSize',14);
        end
    end
end