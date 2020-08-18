%% spm analysis:
spm       = spm1d.stats.ttest2(YA, YB);
spmi      = spm.inference(0.05, 'two_tailed',true, 'interp',true);
% disp(spmi)



%(2) Plot:
% close all
% spmi.plot();
% spmi.plot_threshold_label();
% spmi.plot_p_values();

%% custom spm plot
plot_distribution(time,avethor_flex_mstr_c)
plot_distribution(time,avethor_flex_mstr_y)
ax = gca;
color1 = ax.ColorOrder(1,:);
color2 = ax.ColorOrder(2,:);
color3 = ax.ColorOrder(3,:);
ylims = [min(ylim) min(ylim) max(ylim) max(ylim)];
close all
    
fig1 = figure('color','w');
    subplot(2,1,1)
        box('on')
        sgtitle(["\fontsize{16}Trunk Flexion",...
            '\fontsize{10}\color[rgb]{0 0.447 0.741}college \color{black}| \color[rgb]{0.929 0.694 0.125}youth']...
            ,'fontweight','normal','fontname','times new roman');
        ylabel(['Joint Angle (' char(176) ')'])
        xline(mean(clock3_time_y),'linestyle','--','color',color1)
        xline(mean(clock12_time_y),'linestyle','--','color',color1)
        xline(mean(fc_time_y),'linestyle','--','color',color1)
        xline(mean(br_time_y),'linestyle','--','color',color1)
        xline(mean(clock3_time_c),'linestyle','--','color',color3)
        xline(mean(clock12_time_c),'linestyle','--','color',color3)
        xline(mean(fc_time_c),'linestyle','--','color',color3)
        xline(mean(br_time_c),'linestyle','--','color',color3)
        yline(0);
        % patch(clock3xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(clock12xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(fcxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        % patch(brxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
        plot_distribution(time,avethor_flex_mstr_c);
        plot_distribution(time,avethor_flex_mstr_y);
        text((clock3mean_c+clock3mean_y)/2,max(ylim),'3:00','fontname',...
            'times new roman','fontsize',8,'horizontalalignment','center',...
            'verticalalignment','bottom')
        text((clock12mean_c+clock12mean_y)/2,max(ylim),'12:00','fontname',...
            'times new roman','fontsize',8,'horizontalalignment','center',...
            'verticalalignment','bottom')
        text((fcmean_c+fcmean_y)/2,max(ylim),'FC','fontname','times new roman',...
            'fontsize',8,'horizontalalignment','center','verticalalignment',...
            'bottom')
        text((brmean_c+brmean_y)/2,max(ylim),'BR','fontname','times new roman',...
            'fontsize',8,'horizontalalignment','center','verticalalignment',...
            'bottom')
        set(gca,'fontname','times new roman','fontsize',12);
        xticks([0 .25 .50 .75 1]);
        xticklabels(["" "25" "50" "75" ""]);
    subplot(2,1,2)
        spmi.plot();
        spmi.plot_threshold_label();
        spmi.plot_p_values();
        set(gca,'fontname','times new roman','fontsize',12);
        box('on');
        xlabel('Time (%)')
        xticks([0 25 50 75 100]);
        xticklabels(["" "25" "50" "75" ""]);