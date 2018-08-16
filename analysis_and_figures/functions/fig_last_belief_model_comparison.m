function fig_last_belief_model_comparison(datadir,modeldir)

%%
% datadir = '../../data/data_CCONF.mat';
uni_sujs = 1:3;
[d,m_like] = get_data(uni_sujs,datadir,modeldir,'use_like_nonparW_shrink');
[d,m_ignore] = get_data(uni_sujs,datadir,modeldir,'ignore_conf_nonparW_shrink');
[d,m_conf] = get_data(uni_sujs,datadir,modeldir,'use_conf_nonparW_shrink');

%%
I = d.trnum_back_eff==1;
dife_like = bsxfun(@minus,d.belief(I), m_like.belief_model(I,:));
dife_conf = bsxfun(@minus,d.belief(I), m_conf.belief_model(I,:));
dife_ignore = bsxfun(@minus,d.belief(I), m_ignore.belief_model(I,:));

a = mean(dife_like.^2);
b = mean(dife_conf.^2);
c = mean(dife_ignore.^2);

% [h,pp] = ttest(a,b);
% [h,pp] = ttest(a,c);

%%
[newScat, newHist] = fn_hist_plot(a',b',[0.02,0.045],'MSE, Bayesian model','MSE, Choice-conf model',120,12,0);
pos = get(newHist,'position');
set(newHist,'position',[0.04,0.04,pos(3),pos(4)])

[newScat, newHist] = fn_hist_plot(a',c',[0.02,0.065],'MSE, Bayesian model','MSE, Choice-only model',120,12,0);
pos = get(newHist,'position');
set(newHist,'position',[0.07,0.07,pos(3),pos(4)])

end

