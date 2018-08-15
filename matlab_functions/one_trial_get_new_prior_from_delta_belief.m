function [prior_K,like] = one_trial_get_new_prior_from_delta_belief(K,prior_K_prev,choice,delta_belief)
% K: vector of biases: eg: [0:0.2:1]
% prior: probability of each bias
% choice: 1 right, 0 left
% conf: 1 sure right, 0.5 guessing

% THIS IS CONCEPTUALLY WRONG. DELTA BELIEF MAY NOT CHANGE WHEN PRIOR_K
% DOES, FOR INSTANCE WHEN THE Ss IS ALREADY CERTAIN IT IS BIASED


if min(choice)<0 || delta_belief>1 || delta_belief<-1
    error('problem in input')
end


pRightBiasPrev = marginalize_bias(prior_K_prev,K);

f = @(theta) (calc_diff(theta, delta_belief, choice,K, prior_K_prev, pRightBiasPrev));
like = fminbnd(f, 0, 1);

[err,prior_K] = calc_diff(like, delta_belief, choice,K, prior_K_prev, pRightBiasPrev);

end


function [err,pk] = calc_diff(theta, delta_belief, choice,K, prior_K_prev, pRightBiasPrev)

like = theta(1);

p_right = choice.*like + (1-choice).*(1-like);%evidence
pev_k = K.*p_right+(1-K).*(1-p_right);
pk_ev = pev_k.*prior_K_prev;
pk_ev = pk_ev/sum(pk_ev);
pk = pk_ev;

pRightBias = marginalize_bias(pk,K);

delta_b_pred = pRightBias - pRightBiasPrev;

err = (delta_b_pred-delta_belief).^2;

end

function pRightBias = marginalize_bias(prior_K,K)
pRightBias = sum(prior_K(K>0.5));
if sum(K==0.5)>0
    pRightBias = pRightBias  + sum(prior_K(K==0.5))/2;
end
end
