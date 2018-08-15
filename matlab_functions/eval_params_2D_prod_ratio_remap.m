function [logl, model, probMapUnbiased, P] = eval_params_2D_prod_ratio_remap(theta,D,seed,pars)
% function [logl, model, probMapUnbiased, P] = eval_params_2D_prod_ratio_remap(theta,D,seed,pars)

% wtype:
%1: 'use_like_nonparW_shrink'
%2: 'use_conf_nonparW_shrink'
%3: 'ignore_conf_nonparW_shrink'

kappa   = theta(1);
B0      = theta(2);
coh0    = theta(3);

w2      = theta(4);
w3      = theta(5);
shrink  = theta(6);
phi     = theta(7);
K = [0.6-(0.6-0.5)*shrink,...
    0.8-(0.8-0.5)*shrink,...
    1-(1-0.5)*shrink];
K = sort([1-K,K]);
prior_K_ini = [w3,w2,1,1,w2,w3];
prior_K_ini = prior_K_ini/sum(prior_K_ini);

bupdate = [theta(9:13)];

fitconf_flag = pars.fitconf_flag;

if isfield(pars,'is_for_eval') && pars.is_for_eval==1
    save_extra = 1;
else
    save_extra = 0;
end

%% scale confidence to 0.5-1 range
use_remap = 1;
if use_remap
    
    nsims = 1;
    [~,model_allsim] = simulate_process_2D_prod_ratio(theta,D,K,prior_K_ini,nsims,seed+1);
    
    conf = sort_two_vec_reassign(D.confidence,model_allsim.conf(:,1));

    % save
    model.conf_remapped = conf;
    
else
    conf = clip(D.confidence,-1,1);
    conf = (1+conf)/2;
end

%% dtb, unbiased case

theta_unbiased = [kappa,B0,coh0,0];
Ps = wrapper_dtb_noprior(theta_unbiased,D.dotdur,D.coh, 1 ); %overkill
P  = wrapper_dtb_noprior(theta_unbiased,D.dotdur,D.coh, 0 ); 
signed = sign(Ps.drift)';


%% confidence map unbiased
probMapUnbiased = beliefmap_1d_sdrift(Ps.notabs.pdf,Ps.up.pdf_t,Ps.lo.pdf_t,Ps.y, signed ,ones(size(Ps.drift)));

% add a field to probMap
pr = probMapUnbiased.pcorr_nobound;
I = probMapUnbiased.best_choice_nobound==0;
pr(I) = 1 - pr(I);
probMapUnbiased.pcorr_given_right_nobound = pr;

probMapUnbiased.pcorr_given_right_bound_up = probMapUnbiased.pcorr_bound_up;
probMapUnbiased.pcorr_given_right_bound_lo = 1-probMapUnbiased.pcorr_bound_up;


%%
[~,~,icoh] = nanunique(D.coh);
rng(seed,'twister');

%% simulate

ntr = length(D.choice);
model.p_choice      = nan(ntr,1);
model.p_for_belief_update = nan(ntr,1);
model.p_choice_conf       = nan(ntr,1);
model.e_conf              = nan(ntr,1);
model.prior_rightwards    = nan(ntr,1);
model.right_choice_posterior_unbiased = nan(ntr,1);

model.K = K;
model.prior_K_ini = prior_K_ini;
model.prior_K    = nan(ntr,length(K));
model.posterior_K    = nan(ntr,length(K));


for i = 1:ntr
    
    %first trial of block
    if i==1 || (D.session(i) ~= D.session(i-1)) 
        prior_K = prior_K_ini;
    end
    
    % prob of rightward
    if fitconf_flag
        if D.choice(i)==1
            if D.high_confidence(i)==1
                choice_posterior_range = [phi,1];
            else
                choice_posterior_range = [0.5,phi];
            end
        else
            if D.high_confidence(i)==1
                choice_posterior_range = [0,1-phi];
            else
                choice_posterior_range = [1-phi,0.5];
            end
        end
    else
        if D.choice(i)==1
            choice_posterior_range = [0.5,1];
        else
            choice_posterior_range = [0,0.5];
        end
    end
    
    calc_expected_conf_flag = false;
    [prob_choice_belief_in_range,for_sampling] = p_choice_conf_in_range(icoh(i),D.dotdur(i),...
        choice_posterior_range,probMapUnbiased,P,K,prior_K, calc_expected_conf_flag);

    
    if save_extra
        % just for the evaluation of the best one, not while fitting
        
        choice_range = [0.5,1]*D.choice(i) + [0,0.5]*(1-D.choice(i));
        
        [model.p_choice(i),~,model.e_conf(i)] = p_choice_conf_in_range(icoh(i),D.dotdur(i),...
            choice_range,probMapUnbiased,P,K,prior_K, true );
        
        model.prior_K(i,:) = prior_K;
    end
    
    % update prior for next trial
    prior_rightwards = prior_K*K(:);
    pnext_rel_choice = D.choice(i)*prior_rightwards + (1-D.choice(i))*(1-prior_rightwards);
    yy = bupdate * [1,conf(i),pnext_rel_choice, conf(i)*pnext_rel_choice, conf(i)/pnext_rel_choice]';
    p_for_belief_update = 1./(1+exp(-yy));
    
    % now update belief:
    [prior_K,~] = one_trial_bayesian_update(K,prior_K,D.choice(i),p_for_belief_update);    
    
    
    model.p_for_belief_update(i) = p_for_belief_update;
    model.p_choice_conf(i) = prob_choice_belief_in_range;
    
    model.prior_rightwards(i) = prior_rightwards;
    
    if save_extra
        model.posterior_K(i,:) = prior_K;
    end
end



%% eval like
pPred = model.p_choice_conf;

pPred(pPred<=0) = eps;

logl = -sum(log(pPred));

%% print
fprintf('err=%.3f, kappa=%.2f, B0=%.2f,  coh0=%.2f, w2=%.2f, w3=%.2f, shrink=%.2f, phi=%.2f, b1=%.2f, b2=%.2f, b3=%.2f, b4=%.2f ,b5=%.2f \n',...
    logl , theta(1), theta(2), theta(3), theta(4), theta(5), theta(6), theta(7), theta(9), theta(10), theta(11), theta(12),theta(13));


