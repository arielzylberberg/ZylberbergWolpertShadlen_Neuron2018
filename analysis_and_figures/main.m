% restoredefaultpath % to check dependencies; remove later
clear
addpath(genpath('../matlab_functions/'));
addpath('./functions');

%%

data_path_to_file = '../data/data_CCONF.mat';
modeldir = '../models/fit_choice_and_conf_3_main_models/';

%% figure 2 
p = figure2();

%% belief pred, model comparison
p = evol_belief_3models(data_path_to_file,modeldir);

%% belief pred last trial, model comparison
fig_last_belief_model_comparison(data_path_to_file,modeldir);
drawnow

%% other figures, prepare
model_type = 1;
uni_sujs = 1:3;
switch model_type
    case 1
        str = 'use_like_nonparW_shrink';
    case 2
        str = 'use_conf_nonparW_shrink';
    case 3
        str = 'ignore_conf_nonparW_shrink';
end

[d,m] = get_data(uni_sujs,data_path_to_file,modeldir,str);
pplot = plot_functions(d,m);

%% go
pplot.dbelief_vs_belief();% moving av - in paper submitted to neuron
pplot.dbelief_vs_belief_equispaced();%new for rev
% pplot.dbelief_vs_belief_unsplit_by_choice_split_prev_coh();
%pplot.choice_vs_belief(); %for reviewer
%pplot.choice_vs_belief_equispaced; %for reviewer
pplot.conf_unsigned_per_suj();
pplot.conf_continuous_unsigned_per_suj();
pplot.acc_unsigned();
pplot.plot_weights(3);
%pplot.influence_bias_0coh_option2(); % Figure 8

pplot.confidence_transformation(); % supp fig

% after neuron's reviews
% if isfield(pplot.m,'like_from_conf')
%     pplot.dbelief_vs_like_from_model();
% end

pplot.combined_fig4_opt2();
    
flag_average_simulations = 0;
p = pplot.combined_fig6(flag_average_simulations);

pplot.combined_supp_fig3();
