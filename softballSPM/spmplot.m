%% spm analysis
% add spm folder to path
spm       = spm1d.stats.ttest2(avethor_rot_mstr_c, avethor_rot_mstr_y);
spmi      = spm.inference(0.05, 'two_tailed',true, 'interp',true);
disp(spmi)

%% custom spm plot
plot_distribution(time,avethor_rot_mstr_c)
plot_distribution(time,avethor_rot_mstr_y)
ax = gca;
color1 = ax.ColorOrder(1,:);
color2 = ax.ColorOrder(2,:);
color3 = ax.ColorOrder(3,:);
ylims = [min(ylim) min(ylim) max(ylim) max(ylim)];
close all
    
fig1 = figure('color','w');
    subplot(2,1,1)
        box('on')
        sgtitle(["\fontsize{16}Thorax Axial Rotation",...
            '\fontsize{10}\color[rgb]{0 0.447 0.741}college \color{black}| \color[rgb]{0.929 0.694 0.125}youth']...
            ,'fontweight','normal','fontname','times new roman');
        ylabel(['Joint Angle (' char(176) ')'])
        %xline(mean(clock3_time_y),'linestyle','--','color',color1)
        xline(mean(clock12_time_y),'linestyle','--','color',color1)
        xline(mean(fc_time_y),'linestyle','--','color',color1)
        %xline(mean(br_time_y),'linestyle','--','color',color1)
        %xline(mean(clock3_time_c),'linestyle','--','color',color3)
        xline(mean(clock12_time_c),'linestyle','--','color',color3)
        xline(mean(fc_time_c),'linestyle','--','color',color3)
        %xline(mean(br_time_c),'linestyle','--','color',color3)
        yline(0);
        % patch(clock3xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(clock12xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(fcxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(brxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        plot_distribution(time,avethor_rot_mstr_c);
        plot_distribution(time,avethor_rot_mstr_y);
%         text((clock3mean_c+clock3mean_y)/2,max(ylim),'3:00','fontname',...
%             'times new roman','fontsize',8,'horizontalalignment','center',...
%             'verticalalignment','bottom')
        text((clock12mean_c+clock12mean_y)/2,max(ylim),'12:00','fontname',...
            'times new roman','fontsize',8,'horizontalalignment','center',...
            'verticalalignment','bottom')
        text((fcmean_c+fcmean_y)/2,max(ylim),'FC','fontname','times new roman',...
            'fontsize',8,'horizontalalignment','center','verticalalignment',...
            'bottom')
%         text((brmean_c+brmean_y)/2,max(ylim),'BR','fontname','times new roman',...
%             'fontsize',8,'horizontalalignment','center','verticalalignment',...
%             'bottom')
        ax = gca;
        ax.FontName = 'times new roman';
        ax.FontSize = 12;
        ax.XTick = [0 .25 .50 .75 1];
        ax.XTickLabel = ["" "25" "50" "75" ""];
    subplot(2,1,2)
        spmi.plot();
        spmi.plot_threshold_label();
        spmi.plot_p_values();
        box('on');
        ax = gca;
        ax.FontName = 'times new roman';
        ax.FontSize = 12;
        ax.XTick = [0 25 50 75 100];
        ax.XTickLabel = ["" "25" "50" "75" ""];
        xlabel("Time (% 3:00 - BR)");