function [prob_choice_conf_in_range, for_sampling, mean_conf_in_range] = p_choice_conf_in_range(icoh,dotdur,...
    choice_posterior_range,probMapUnbiased,P,K,prior_K, calc_expected_conf_flag)



prior_K = prior_K/sum(prior_K); %norm

t = P.t;
dt = t(2)-t(1);
isamp = round(dotdur/dt);

% no bound
pr = probMapUnbiased.pcorr_given_right_nobound(:,isamp);
pcorr_given_right = prior_K*K(:)*pr;
pcorr_given_left = prior_K*(1-K(:))*(1-pr);
pcorr_given_right = pcorr_given_right./(pcorr_given_right + pcorr_given_left);

% bound up
pr = probMapUnbiased.pcorr_given_right_bound_up;
pcorr_bound_up_given_right = prior_K*K(:)*pr;
pcorr_bound_up_given_left = prior_K*(1-K(:))*(1-pr);
pcorr_bound_up_given_right = pcorr_bound_up_given_right./(pcorr_bound_up_given_right + pcorr_bound_up_given_left);

% bound lo
pr = probMapUnbiased.pcorr_given_right_bound_lo;
pcorr_bound_lo_given_right = prior_K*K(:)*pr;
pcorr_bound_lo_given_left = prior_K*(1-K(:))*(1-pr);
pcorr_bound_lo_given_right = pcorr_bound_lo_given_right./(pcorr_bound_lo_given_right + pcorr_bound_lo_given_left);


% indexes
f_choice_notabs = pcorr_given_right>=choice_posterior_range(1) & pcorr_given_right<=choice_posterior_range(2);
f_choice_up = pcorr_bound_up_given_right>=choice_posterior_range(1) & pcorr_bound_up_given_right<=choice_posterior_range(2);
f_choice_lo = pcorr_bound_lo_given_right>=choice_posterior_range(1) & pcorr_bound_lo_given_right<=choice_posterior_range(2);

f_notabs = f_choice_notabs;
f_up = f_choice_up;
f_lo = f_choice_lo;

prob_choice_conf_in_range = sum(P.notabs.pdf(icoh,f_notabs,isamp)) + sum(P.up.pdf_t(icoh,f_up(1:isamp))) + sum(P.lo.pdf_t(icoh,f_lo(1:isamp)));


% for sampling
for_sampling.lo_pdf_t = P.lo.pdf_t(icoh,1:isamp).*f_lo(1:isamp);
for_sampling.up_pdf_t = P.up.pdf_t(icoh,1:isamp).*f_up(1:isamp);
for_sampling.notabs_pdf = P.notabs.pdf(icoh,:,isamp).*f_notabs';

if calc_expected_conf_flag
    pcorr = max(pcorr_given_right,1-pcorr_given_right);
    pcorr_bound_up = max(pcorr_bound_up_given_right,1-pcorr_bound_up_given_right);
    pcorr_bound_lo = max(pcorr_bound_lo_given_right,1-pcorr_bound_lo_given_right);
    
    aux = P.notabs.pdf(icoh,f_notabs,isamp)*pcorr(f_notabs) + ...
        P.up.pdf_t(icoh,f_up(1:isamp))*pcorr_bound_up(f_up(1:isamp))' + ...
        P.lo.pdf_t(icoh,f_lo(1:isamp))*pcorr_bound_lo(f_lo(1:isamp))';
    
    mean_conf_in_range = aux/prob_choice_conf_in_range;
    
else
    mean_conf_in_range = [];
    
end

end
