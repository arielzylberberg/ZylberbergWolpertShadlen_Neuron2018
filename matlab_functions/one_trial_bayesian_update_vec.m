function [Ph,pRightBias] = one_trial_bayesian_update_vec(K,prior,v_choice,v_conf)
% K: vector of biases: eg: [0:0.2:1]
% prior: probability of each bias
% choice: 1 right, 0 left
% conf: 1 sure right, 0.5 guessing

if min(v_choice)<0 || max(v_conf)>1 || min(v_conf)<0
    error('problem in input')
end

% N = length(choice);

%start with uniform prior over bias
ph = prior;

p_right = v_choice.*v_conf + (1-v_choice).*(1-v_conf);%evidence

%p(H/ev) = p(ev/H)*p(H)/p(ev)

% for i=1:N
ev = p_right;
pev_h = K(:)*ev(:)'+(1-K(:))*(1-ev(:)');
% ph_ev = pev_h.*ph;
ph_ev = bsxfun(@times,pev_h,ph(:));
ph_ev = bsxfun(@times,ph_ev,1./sum(ph_ev));
Ph = ph_ev;
%update prior
%     ph = ph_ev;
% end

pRightBias = sum(Ph(K>0.5,:));
if sum(K==0.5)>0
    pRightBias = pRightBias  + sum(Ph(K==0.5,:))/2;
end

end

