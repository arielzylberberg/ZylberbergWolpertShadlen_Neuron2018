function [logl,model_allsim, P, suave] = simulate_process2(theta,D,K,prior_K_ini,nsims,seed,smooth_flag,wtype)
% function [logl,model_allsim, P, suave] = simulate_process2(theta,D,K,prior_K_ini,nsims,seed,smooth_flag,wtype)

% wtype:
%1: 'use_like_nonparW_shrink'
%2: 'use_conf_nonparW_shrink'
%3: 'ignore_conf_nonparW_shrink'

kappa   = theta(1);
B0      = theta(2);
coh0    = theta(3);
w2      = theta(4);
w3      = theta(5);
shrink = theta(6);
phi = theta(7);

prior_K_ini = prior_K_ini/sum(prior_K_ini);

if isequal(wtype,'ignore_conf_nonparW_shrink') || wtype == 3
    fix_like = theta(8);
end

%% dtb, unbiased case

theta_unbiased = [kappa,B0,coh0,0];
Ps = wrapper_dtb_noprior(theta_unbiased,D.dotdur,D.coh, 1 ); %overkill
P = wrapper_dtb_noprior(theta_unbiased,D.dotdur,D.coh, 0 ); 

%% sample dv and t for all trials - takes forever

rng(seed,'twister');

[ucoh,~,icoh] = unique(D.coh);
ntr = length(D.dotdur);
idv      = nan(ntr,nsims);
samp_num = nan(ntr,nsims);
for i=1:nsims
    [idv(:,i),samp_num(:,i)] =  sample_dv_t_faster(P,D.dotdur,icoh);
end

%% for each sample, get the probability that a rightward response is 
% correct in the unbiased condition
right_conf_nobias = calc_pright_resp_mbote(Ps.y(idv),Ps.t(samp_num),Ps.drift);
if isrow(right_conf_nobias)
    right_conf_nobias = right_conf_nobias';
end

%% simulate many runs per block

model_allsim.choice = nan(ntr,nsims);
model_allsim.like   = nan(ntr,nsims);
model_allsim.conf   = nan(ntr,nsims);
model_allsim.belief = nan(ntr,nsims);
nk = length(K);
model_allsim.posterior_belief = nan(ntr,nsims,nk);
model_allsim.prior_belief = nan(ntr,nsims,nk);

for j=1:nsims
    for i = 1:ntr
        if i==1 || (D.session(i) ~= D.session(i-1)) 
            prior_K = prior_K_ini;
        end
        
        pR_g_et = (prior_K * K') * right_conf_nobias(i,j);
        pL_g_et = (prior_K * (1-K)') * (1-right_conf_nobias(i,j));
        pR_g_et = pR_g_et/(pR_g_et+pL_g_et);
        pL_g_et = 1 - pR_g_et;
        
        if pR_g_et>0.5
            m_choice = 1;
            m_conf = pR_g_et;
        else
            m_choice = 0;
            m_conf = 1-pR_g_et;
        end
        
        if wtype == 2 || isequal(wtype,'use_conf_nonparW_shrink')
            m_like = m_conf;
        elseif wtype == 3 || isequal(wtype,'ignore_conf_nonparW_shrink')
            m_like = fix_like;
        else
            m_like = m_choice*right_conf_nobias(i,j) + (1-m_choice)*(1-right_conf_nobias(i,j));
        end
        
        prior_K_prev = prior_K;
        [prior_K,~] = one_trial_bayesian_update(K,prior_K, m_choice ,m_like);
        
        % belief
        PRIGHTBIAS = sum(prior_K(K>0.5));
        if any(K==0.5)
            PRIGHTBIAS = PRIGHTBIAS + sum(prior_K(K==0.5))/2;
        end
        
        %
        model_allsim.choice(i,j) = m_choice;
        model_allsim.like(i,j)   = m_choice*right_conf_nobias(i,j) + (1-m_choice)*(1-right_conf_nobias(i,j));
        model_allsim.conf(i,j)   = m_conf;
        model_allsim.belief(i,j) = PRIGHTBIAS;
        model_allsim.posterior_belief(i,j,:) = prior_K; %save the full posterior
        model_allsim.prior_belief(i,j,:) = prior_K_prev;%overkill
        
    end
end


%% smooth lines

logl = nan;

% %% 
if smooth_flag

    ucohs = linspace(eps,0.6,60);
    ucohs = sort([-ucohs,ucohs]);
    
    
    theta_unbiased = [kappa,B0,coh0,0];
    Pf = wrapper_dtb_noprior(theta_unbiased,D.dotdur,ucohs, 0 );
    
    signed = sign(Ps.drift)';
    probMapUnbiased = beliefmap_1d_sdrift(Ps.notabs.pdf,Ps.up.pdf_t,Ps.lo.pdf_t,Ps.y, signed ,ones(size(Ps.drift)));
    
    % add a field to probMap
    pr = probMapUnbiased.pcorr_nobound;
    I = probMapUnbiased.best_choice_nobound==0;
    pr(I) = 1 - pr(I);
    probMapUnbiased.pcorr_given_right_nobound = pr;
    
    probMapUnbiased.pcorr_given_right_bound_up = probMapUnbiased.pcorr_bound_up;
    probMapUnbiased.pcorr_given_right_bound_lo = 1-probMapUnbiased.pcorr_bound_up;
    
    p_right = nan(ntr,length(ucohs));
    p_right_high = nan(ntr,length(ucohs));
    p_left_high = nan(ntr,length(ucohs));
    econf_right = nan(ntr,length(ucohs));
    econf_left = nan(ntr,length(ucohs));
    for j=1:ntr
        prior_K = squeeze(model_allsim.prior_belief(j,1,:))';
        for i=1:length(ucohs)

            p_right(j,i) = p_choice_conf_in_range(i,D.dotdur(j),...
                [0.5,1],probMapUnbiased,Pf,K,prior_K, 0);

            p_right_high(j,i) = p_choice_conf_in_range(i,D.dotdur(j),...
                [phi,1],probMapUnbiased,Pf,K,prior_K, 0);

            p_left_high(j,i) = p_choice_conf_in_range(i,D.dotdur(j),...
                [0,1-phi],probMapUnbiased,Pf,K,prior_K, 0);
            
            [~,~,econf_right(j,i)] = p_choice_conf_in_range(i,D.dotdur(j),...
                [0.5,1],probMapUnbiased,Pf,K,prior_K, 1);
            
            [~,~,econf_left(j,i)] = p_choice_conf_in_range(i,D.dotdur(j),...
                [0,0.5],probMapUnbiased,Pf,K,prior_K, 1);
        end
    end
    
    % With more durations, representative sampling
    sdotdur = sort(D.dotdur);
    I = ceil(linspace(0,1,20)*length(sdotdur));
    I(1) = 1;
    udurs = sdotdur(I);
    
    p_right_d = nan(ntr,length(ucoh),length(udurs));
    p_right_high_d = nan(ntr,length(ucoh),length(udurs));
    p_left_high_d = nan(ntr,length(ucoh),length(udurs));
    
    econf_right_d = nan(ntr,length(ucoh),length(udurs));
    econf_left_d = nan(ntr,length(ucoh),length(udurs));
    for j=1:ntr
        prior_K = squeeze(model_allsim.prior_belief(j,1,:))';
        for i=1:length(Ps.drift) % or P?
            for k=1:length(udurs)
                p_right_d(j,i,k) = p_choice_conf_in_range(i,udurs(k),...
                    [0.5,1],probMapUnbiased,Ps,K,prior_K, 0);
                
                p_right_high_d(j,i,k) = p_choice_conf_in_range(i,udurs(k),...
                    [phi,1],probMapUnbiased,Ps,K,prior_K, 0);
                
                p_left_high_d(j,i,k) = p_choice_conf_in_range(i,udurs(k),...
                    [0,1-phi],probMapUnbiased,Ps,K,prior_K, 0);
                
                [~,~,econf_right_d(j,i,k)] = p_choice_conf_in_range(i,udurs(k),...
                    [0.5,1],probMapUnbiased,Ps,K,prior_K, 1);
                [~,~,econf_left_d(j,i,k)] = p_choice_conf_in_range(i,udurs(k),...
                    [0.0,0.5],probMapUnbiased,Ps,K,prior_K, 1);
            end
        end
    end
    
    suave = struct('p_right',p_right,'p_right_high',p_right_high,'p_left_high',p_left_high,'ucohs',ucohs,...
        'p_right_d',p_right_d,'p_right_high_d',p_right_high_d,'p_left_high_d',p_left_high_d,'udurs',udurs,...
        'econf_right',econf_right,'econf_left',econf_left,'econf_right_d',econf_right_d,'econf_left_d',econf_left_d);
    
    figure();
    choice_model = nanmean(model_allsim.choice,2);
    u = unique(D.block_prior);
    colores = colors_az(length(u));
    for i=1:length(u)
        I = D.block_prior==u(i);
        plot(ucohs,nanmean(suave.p_right(I,:),1),'color',colores(i,:));
        hold all
        [tt,yy] = curva_media(choice_model,D.coh,I,0);% sanity check
        plot(tt,yy,'color',colores(i,:));
        [tt,yy,ss] = curva_media(D.choice,D.coh,I,0);
        terrorbar(tt,yy,ss,'color',colores(i,:),'marker','o','LineStyle','none','markerfacecolor',colores(i,:),'markersize',6);
    end
else
    suave = [];
end


%% print
fprintf('err=%.3f, kappa=%.2f, B0=%.2f,  coh0=%.2f, w2=%.2f, w3=%.2f, shrink=%.2f, phi=%.2f\n',...
    logl , theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7));


