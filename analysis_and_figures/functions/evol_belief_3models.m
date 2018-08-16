function p = evol_belief_3models(datadir,modeldir)

uni_sujs = 1:3;
folders = {'use_like_nonparW_shrink','use_conf_nonparW_shrink','ignore_conf_nonparW_shrink'};
str = {'Bayesian model','Choice-Confidence model','Choice-only model'};

p = publish_plot(3,3);

for i=1:length(folders)

    [d,m] = get_data(uni_sujs,datadir,modeldir,folders{i});

    pplot = plot_functions(d,m);
    
    p_evol = pplot.belief_vs_trnum2();
    p.copy_from_ax(p_evol.h_ax,i);
    ht = title(str{i});
    set(ht,'fontweight','normal');
    
    pp = pplot.combined_fig6(0);
    p.copy_from_ax(pp.h_ax(8),i + 3);
    
    p.copy_from_ax(pp.h_ax(9),i + 6);
    
    close(p_evol.h_fig)
    close(pp.h_fig);
    
end

[idx,I1,I2] = p.center_plots();
for i=1:length(idx)
    p.current_ax(idx(i));
    ylabel('');
    xlabel('');
    set(gca,'yticklabel','','xticklabel','');
end

p.current_ax(8);ylabel('');
p.current_ax(9);ylabel('');

set(p.h_ax(7:9),'xtick',[0:0.5:1],'xticklabel',0:0.5:1)



for i=4:6
    ch = findall(get(p.h_ax(i),'children'),'Type','Line');
    set(ch,'markerfacecolor','k','markeredgecolor','w','markersize',6);
end


% p.unlabel_center_plots();
% set(gcf,'Position',[401   93  599  705])
set(gcf,'Position',[482  384  599  554])
p.format('FontSize',11);

drawnow

end

