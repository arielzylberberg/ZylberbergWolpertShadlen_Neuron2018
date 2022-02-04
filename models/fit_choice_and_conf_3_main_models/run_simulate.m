% restoredefaultpath % to check dependencies; remove later
clear
addpath(genpath('../../matlab_functions/'));

%% load data
load('../../data/data_CCONF.mat');
nsuj = length(dat);

%% combinations of subjects and fit types
combs = combvec([1:nsuj],[1,2,3]); %subject,mtype

%%

parfor i=1:size(combs,1)
% for i=1:size(combs,1)
    
    suj = combs(i,1);
    wtype = combs(i,2);
    
    pars = struct();
    switch wtype
        case 1
            pars.wtype = 'use_like_nonparW_shrink';
        case 2
            pars.wtype = 'use_conf_nonparW_shrink';
        case 3
            pars.wtype = 'ignore_conf_nonparW_shrink';
    end
    
    datadir = pwd;
    filename = fullfile(datadir,pars.wtype,['suj_',num2str(suj)]);
    mfits = load(filename);
    
    D = dat(suj);
    
    %% make K and prior_K_ini
    w2      = mfits.theta(4);
    w3      = mfits.theta(5);
    shrink  = mfits.theta(6);
    
    K = [0.6-(0.6-0.5)*shrink,...
        0.8-(0.8-0.5)*shrink,...
        1-(1-0.5)*shrink];
    K = sort([1-K,K]);
    prior_K_ini = [w3,w2,1,1,w2,w3];
    
    
    %% perform simulations with opt parameter
    % nsims = 200; % THIS IS 200 IN THE PUBLISHED_PAPER
    nsims = 20;
    seed  = 12345;
    smooth_flag = 1;
    [~,model_allsim, P, smooth] = simulate_process2(mfits.theta,D,K,prior_K_ini,nsims,seed ,smooth_flag,wtype);
        
    %% reformatting
    f = fields(model_allsim);
    model = struct();
    for j=1:length(f)
        model.(f{j}) = nanmean(model_allsim.(f{j}),2);
    end
    for j=1:length(f)
        model.(f{j}) = to_mat(model.(f{j}),dat(suj).session);
    end
    
    struct_to_save = struct('model',model,'model_allsim',model_allsim,'smooth',smooth);
    save_parallel(filename,struct_to_save, 1 );
    
    
end

