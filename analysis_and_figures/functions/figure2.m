function pplot = figure2()


%%

uni_sujs = 1:3;

datadir = '../data/data_CCONF.mat';

[d,m] = get_data(uni_sujs,datadir);


%%

struct2vars(d);

nSujs = length(uni_sujs);

uni_block_bias = unique(block_prior_rel);

pplot = publish_plot(2,2);
set(gcf,'Position',[383  241  679  557])
pplot.new_axes('Position',[0.37,0.66,0.15,0.15]);
pplot.displace_ax(2,0.06,1)
pplot.displace_ax(4,0.06,1)

% set(gcf,'Position',[469  655  431  345])

markersize = 7;
markersize_logplot = 9;
edgecolor = 'none';

%% do regression

filt = true(size(choice));

ntr = sum(filt);
depvar = choice(filt);
dummy_prior = adummyvar(block_prior);
dummy_group = adummyvar(group);
indepvar = {'coh_prior',bsxfun(@times,coh(filt),dummy_prior(filt,:)),...
    'dot_dur_coh',bsxfun(@times,dotdur(filt),coh(filt))...
    'prior',dummy_prior(filt,:),...
    'group',dummy_group(filt,1:end-1)};
testSignificance.vars = [1,2,3,4];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

yhat = glmval(beta,x,'logit','constant','off');

% stats
fprintf(' Coh stats: %2.7f\n',LRT(1).p)
fprintf(' Prior stats: %2.7f\n',LRT(3).p)


%% eval the regression model w/ a finer set of coherences, for plotting
coh_fine = linspace(-0.6,0.6,101);

unib = unique(block_prior);
clear yhat_fine
for j=1:length(unib)
    dummy_prior_aux = zeros(ntr,6);
    dummy_prior_aux(:,j) = 1;
    for i=1:length(coh_fine)
        cc = coh_fine(i);
        aux = {'coh_prior',bsxfun(@times,cc*ones(ntr,1),dummy_prior_aux),...
            'dot_dur_coh',bsxfun(@times,dotdur(filt),cc*ones(ntr,1))...
            'prior',dummy_prior_aux,...
            'group',dummy_group(filt,1:end-1)};
        xx = cat(2,aux{2:2:end});
        yhat_fine(:,i,j) = glmval(beta,xx,'logit','constant','off');
    end
end



%% choice vs coh
pplot.next();
[p,X,Y,S] = curva_media_wrapper(choice(filt),coh(filt),block_prior(filt),0);

colores_prior = cbrewer('qual','Paired',6);
colores_prior = colores_prior([1,3,5,6,4,2],:);

for i=1:6
    if unib(i)==0
        I = coh_fine<=0;
    elseif unib(i)==1
        I = coh_fine>=0;
    else
        I = true(size(coh_fine));
    end
    hp(i) = plot(coh_fine(I),nanmean(yhat_fine(:,I,i)),'color',colores_prior(i,:),'linestyle','-','linewidth',1);
    hold all
end

for i=1:length(X)
    terrorbar(X{i},Y{i},S{i},'color',colores_prior(i,:),'marker','o','markersize',markersize,'linestyle','none',...
        'markerfacecolor',colores_prior(i,:),'markeredgecolor',edgecolor);
    hold all
end



[hl(1),icons]  = legend(hp(end:-1:1),{'1','0.8','0.6','0.4','0.2','0'});
ht(1) = get(hl(1),'title');
% set(ht(1),'String','Base rate');
set(hl(1),'location','northwest');
for i=1:length(icons)
    if isequal(icons(i).Type,'line') && length(get(icons(i),'XData'))==2
        set(icons(i),'XData',[0.2,0.5]);
    end
end



set(gca,'xlim',[-0.6,0.6]);
xlabel('Motion coherence')
ylabel('Proportion rightward choices')

%% regression inset
pplot.current_ax(5);
color = 'k';

terrorbar(unib,beta(idx.prior),...
    stats.se(idx.prior),'marker','none','LineStyle','none','color',color);
hold all
lsline

msize = 6;
plot(unib,beta(idx.prior),'linestyle','none','markerfacecolor',color,'markeredgecolor',edgecolor,'markersize',msize,'marker','o')

pxlim
axis tight
symmetric_y(gca)
set(gca,'xtick',[0,0.5,1]);
xlabel('Base-rate of block')

ylabel('Leverage on choice')


%% accuracy vs coh

pplot.current_ax(2);

color = movshon_colors(length(uni_block_bias));


% color = colores_prior([4:6],:);

for i=1:length(uni_block_bias)
    J = block_prior_rel==uni_block_bias(i) & filt;
    [~,hplot(i)] = dotsanalysis.plot_log(coh(J),correct(J),coh(J),correct(J),'linestyle','-','color',color(i,:),'marker','o','markersize',markersize_logplot,'break_scaling',0,...
        'markeredgecolor',edgecolor);
    hold all
    
end

xlabel('Motion strength (%coh)')
ylabel('Proportion correct')

[hl(2),icons]  = legend(hplot(end:-1:1),{'1, 0','0.8, 0.2','0.6, 0.4'});
% ht(2) = get(hl(2),'title');
% set(ht(2),'String','Bias strength');
set(hl(2),'location','southeast');

for i=1:length(icons)
    if isequal(icons(i).Type,'line') && length(get(icons(i),'XData'))==2
        set(icons(i),'XData',[0.1,0.3]);
    end
end

BreakXAxis(gca,0.18,'break_scaling',3);


%% conf vs coh

pplot.current_ax(3);

confidence = 0.5 + conf/2;

[tt,xx,ss] = curva_media_hierarch(confidence,coh,bl,correct==1,0);

color = movshon_colors(5);
for i=1:5
    h(i) = plot(tt,xx(:,i),'color',color(i,:));
    hold all
end
for i=1:5
    terrorbar(tt,xx(:,i),ss(:,i),'color',color(i,:),'marker','o','markersize',markersize,'linestyle','none','markerfacecolor',color(i,:),...
        'markeredgecolor',edgecolor);
end

include_errors_conf = 0;
if include_errors_conf
    [tte,xxe,sse] = curva_media_hierarch(confidence,coh,bl,trnum_eff>=min_tr_subset2 & correct==0,0);
    for i=1:5
        plot(tte,xxe(:,i),'color',color(i,:),'LineStyle','--');
        terrorbar(tte,xxe(:,i),sse(:,i),'color',color(i,:),'marker','o','markersize',markersize,'linestyle','none','markerfacecolor','none',...
            'markeredgecolor',color(i,:));
    end
end


ubl = unique(bl);

[hl(3),icons] = legend(h(end:-1:1),{'1','0.8','0.6','0.4','0.2'});

set(hl(3),'location','southeast');

for i=1:length(icons)
    if isequal(icons(i).Type,'line') && length(get(icons(i),'XData'))==2
        set(icons(i),'XData',[0.2,0.5]);
    end
end


% legend('early trials','late trials');
ylim([0.5,1])

xlabel('Motion coherence')
ylabel('Confidence')

xlim([-0.6,0.6]);



%% belief vs trial-number

pplot.current_ax(4);


for i=1:length(unib)
    I = block_prior==unib(i);
    [tt,xx,ss] = curva_media( belief ,trnum_eff,I,0);
    %[~,hnice(i)] = niceBars2(tt, xx, ss, colores_prior(i,:), .7);
    [~,hnice(i)] = niceBars2(tt, xx, ss, colores_prior(i,:), 1);
    hold all
end
xlim([0,42])
xlabel('Trial number')
ylabel('Belief')


%% pretify the plot

set(pplot.h_ax,'color','none');

set(pplot.h_ax([1,3]),'xtick',[-0.5:0.25:0.5],'xticklabel',[-50:25:50]);

set(pplot.h_ax,'XMinorTick','off','YMinorTick','off');

for i=1:length(pplot.h_ax)
    set(get(pplot.h_ax(5),'children'),'linewidth',1);
end

pplot.format('FontSize',15);
set(pplot.h_ax(5),'FontSize',12);

set(ht,'FontSize',12,'FontWeight','normal');
set(hl,'FontSize',12,'box','off');

%%
drawnow


end


