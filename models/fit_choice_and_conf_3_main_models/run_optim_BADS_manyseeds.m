addpath(genpath('../../matlab_functions/'));

%%

load('../../data/data_CCONF.mat');

nsuj = length(dat);
seed = 22;

%%

% kappa, B0, coh0, w1, w2, shrink, phi, fix_like
tl   = [5  , .3 , 0, 0, 0  , 0, 0.5, 1];
th   = [30 , 5 , 0, 5, 5 , 0, 1, 1];
tg   = [12, 1.5, 0, .5, .5 , 0, 0.75, 1];

Nseeds = 30;
TG = sample_tguess(tl',th',Nseeds,seed + 25235)';

%% combinations of subjects and fit types
combs = combvec([1:nsuj],[1:3],1:Nseeds); %subject,mtype

%%
rng(7891,'twister');

n = size(combs,1);
vtheta = nan(n,length(tl));
vfval = nan(n,1);

to_eval = @(theta,D,pars,wtype) (eval_params_remapping2(theta,D,seed,pars,wtype));

low_prc = 30;
fitconf_flag = 1;

% parfor i=1:size(combs,1)
for i=1:size(combs,1)
    
    suj = combs(i,1);
    wtype = combs(i,2);
    iseed = combs(i,3);
    
    pars = struct('low_prc',low_prc,'fitconf_flag',fitconf_flag);
    
    D = dat(suj);
    D.high_confidence = D.confidence > prctile(D.confidence,pars.low_prc);
    
    
    switch wtype
        case 1
            pars.wtype = 'use_like_nonparW_shrink';
        case 2
            pars.wtype = 'use_conf_nonparW_shrink';
        case 3
            pars.wtype = 'ignore_conf_nonparW_shrink';
    end
    
    %% fit
    fn_fit = @(theta) (to_eval(theta,D,pars,wtype));
    ptl = tl;
    pth = th;
    
    tg = TG(iseed,:);
    
    [vtheta(i,:), vfval(i), exitflag, output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth);
    
end

%% eval best and save
for i = 1:nsuj
    for wtype=1:3
        
        S = find(combs(:,1)==i & combs(:,2)==wtype);
        
        % get best per suj and save
        [fval,I] = min(vfval(S));
        theta = vtheta(S(I),:);
        
        pars = struct('low_prc',low_prc,'fitconf_flag',fitconf_flag,'is_for_eval',1);
        
        suj = i;
        D = dat(suj);
        D.high_confidence = D.confidence > prctile(D.confidence,pars.low_prc);
        
        %% eval
        
        switch wtype
            case 1
                pars.wtype = 'use_like_nonparW_shrink';
            case 2
                pars.wtype = 'use_conf_nonparW_shrink';
            case 3
                pars.wtype = 'ignore_conf_nonparW_shrink';
        end
        
        %eval best
        [~,model_fromeval] = to_eval(theta,D,pars,wtype);
        
        % save
        tosave = struct('theta',theta,'fval',fval,'tl',tl,'th',th,'TG',TG,'vfval',vfval(S),...
            'vtheta',vtheta,'pars',pars,'model_fromeval',model_fromeval);
        
        
        % save
        if ~exist(pars.wtype,'dir')
            mkdir(pars.wtype)
        end
        savefilename = fullfile(pwd,pars.wtype,['suj_',num2str(suj)]);
        save_parallel(savefilename,tosave);
        
    end
end
