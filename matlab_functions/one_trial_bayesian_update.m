function [Ph,pRightBias] = one_trial_bayesian_update(K,prior,choice,conf)
% K: vector of biases: eg: [0:0.2:1]
% prior: probability of each bias
% choice: 1 right, 0 left
% conf: 1 sure right, 0.5 guessing

% sanity check
if min(choice)<0 || max(conf)>1 || min(conf)<0
    error('problem in input')
end


ph = prior;

p_right = choice.*conf + (1-choice).*(1-conf);%evidence

%p(H/ev) = p(ev/H)*p(H)/p(ev)

ev = p_right;
pev_h = K.*ev+(1-K).*(1-ev);
ph_ev = pev_h.*ph;
ph_ev = ph_ev/sum(ph_ev);
Ph = ph_ev;


pRightBias = sum(Ph(K>0.5));
if sum(K==0.5)>0
    pRightBias = pRightBias  + sum(Ph(K==0.5))/2;
end

end

