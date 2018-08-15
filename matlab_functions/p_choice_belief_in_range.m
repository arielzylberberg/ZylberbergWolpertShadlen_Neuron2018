function [prob_choice_belief_in_range, for_sampling, mean_conf_in_range] = p_choice_belief_in_range(icoh,dotdur,...
    choice_posterior_range,belief_posterior_range,probMap,P,K,prior_K,probMapUnbiased, calc_expected_conf_flag)

if nargin<10 || isempty(calc_expected_conf_flag)
    calc_expected_conf_flag = false;
end

prior_K = prior_K/sum(prior_K); %norm

m_like_right_separatrix = nan(1,2);
for i=1:2
    m_like_right_separatrix(i) = find_like_separatrix_given_belief_and_prior(K,prior_K,belief_posterior_range(i));
end

pcorr_given_right = zeros(length(P.y),1);
pcorr_bound_up_given_right = zeros(size(probMap(1).pcorr_bound_up));
pcorr_bound_lo_given_right = zeros(size(probMap(1).pcorr_bound_lo));

t = P.t;
dt = t(2)-t(1);
isamp = round(dotdur/dt);
for i = 1:length(probMap)
    pr = probMap(i).pcorr_given_right_nobound(:,isamp);
    pcorr_given_right = pcorr_given_right + prior_K(i)*pr;%  weighted average
    pcorr_bound_up_given_right = pcorr_bound_up_given_right + prior_K(i)*probMap(i).pcorr_bound_up;
    pcorr_bound_lo_given_right = pcorr_bound_lo_given_right + prior_K(i)*(1-probMap(i).pcorr_bound_lo);
end


%note: there is a numerical issue with the extreme of the interval; I could rewrite to make it 0.5 for each
%category 
f_choice_notabs = pcorr_given_right>=choice_posterior_range(1) & pcorr_given_right<=choice_posterior_range(2);
f_choice_up = pcorr_bound_up_given_right>=choice_posterior_range(1) & pcorr_bound_up_given_right<=choice_posterior_range(2);
f_choice_lo = pcorr_bound_lo_given_right>=choice_posterior_range(1) & pcorr_bound_lo_given_right<=choice_posterior_range(2);


if all(belief_posterior_range==[0,1]) %full range
    f_notabs = f_choice_notabs;
    f_up = f_choice_up;
    f_lo = f_choice_lo;
else
    f_bias_notabs = probMapUnbiased.pright_nobound(:,isamp)>=m_like_right_separatrix(1) & probMapUnbiased.pright_nobound(:,isamp)<=m_like_right_separatrix(2);
    f_bias_up = probMapUnbiased.pright_bound_up>=m_like_right_separatrix(1) & probMapUnbiased.pright_bound_up<=m_like_right_separatrix(2);
    f_bias_lo = probMapUnbiased.pright_bound_lo>=m_like_right_separatrix(1) & probMapUnbiased.pright_bound_lo<=m_like_right_separatrix(2);

    % combined filters
    f_notabs = f_bias_notabs & f_choice_notabs;
    f_up = f_bias_up & f_choice_up;
    f_lo = f_bias_lo & f_choice_lo;
end


prob_choice_belief_in_range = sum(P.notabs.pdf(icoh,f_notabs,isamp)) + sum(P.up.pdf_t(icoh,f_up(1:isamp))) + sum(P.lo.pdf_t(icoh,f_lo(1:isamp)));


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
    
    mean_conf_in_range = aux/prob_choice_belief_in_range;
    
else
    mean_conf_in_range = [];
    
end

end

% if choice_decision==0
%     pcorr_given_choice = 1-pcorr_given_right;
%     pcorr_bound_up_given_choice = 1 - pcorr_bound_up_given_right;
%     pcorr_bound_lo_given_choice = 1 - pcorr_bound_lo_given_right;
% else
%     pcorr_given_choice = pcorr_given_right;
%     pcorr_bound_up_given_choice = pcorr_bound_up_given_right;
%     pcorr_bound_lo_given_choice = pcorr_bound_lo_given_right;
% end


% if choice_belief==1
%     fnotabs = probMapUnbiased.pright_nobound(:,isamp)>m_like_right_separatrix;
%     fup = probMapUnbiased.pright_bound_up>m_like_right_separatrix;
%     flo = probMapUnbiased.pright_bound_lo>m_like_right_separatrix;
% else
%     fnotabs = probMapUnbiased.pright_nobound(:,isamp)<m_like_right_separatrix;
%     fup = probMapUnbiased.pright_bound_up<m_like_right_separatrix;
%     flo = probMapUnbiased.pright_bound_lo<m_like_right_separatrix;
% end

% curtail the distributions
% p_notabs_given_choice_belief = P.notabs.pdf(icoh,:,isamp).*(pcorr_given_choice>0.5)'.*(fnotabs');
% p_up_given_choice_belief = P.up.pdf_t(icoh,:).*(pcorr_bound_up_given_choice>0.5).*(fup);
% p_lo_given_choice_belief = P.lo.pdf_t(icoh,:).*(pcorr_bound_lo_given_choice>0.5).*(flo);
% p_up_given_choice_belief(isnan(p_up_given_choice_belief)) = 0;
% p_lo_given_choice_belief(isnan(p_lo_given_choice_belief)) = 0;


% for each time
% Pdummy.up.pdf_t = p_up_given_choice_belief;
% Pdummy.up.cdf_t = nancumsum(p_up_given_choice_belief,2);
% Pdummy.lo.pdf_t = p_lo_given_choice_belief;
% Pdummy.lo.cdf_t = nancumsum(p_lo_given_choice_belief,2);
% Pdummy.notabs.pdf = p_notabs_given_choice_belief;
% Pdummy.t = P.t;
% Pdummy.y = P.y;
% Pdummy.Bup = P.Bup;
% Pdummy.Blo = P.Blo;
% 
% % dotdur = ones(1000,1)*0.6;
% % icoh = ones(1000,1)*7;
% 
% % sanity check + hack: if the choice is impossible
% impossible = (sum(p_up_given_choice_belief)+sum(p_lo_given_choice_belief)+sum(p_notabs_given_choice_belief(:)))==0;
% if impossible
%     idv = ceil(length(P.y)/2);
%     code = 4;%hack code
%     m_like = 0.5;
%     m_conf = 0;
%     m_choice_unbiased = 1-choice_decision;
% 	return
% end
%   
% 
% [idv,isamp,code] =  sample_one(Pdummy, dotdur);
% 
% %% confidence
% switch code
%     case 1 %bup
%         m_like = probMapUnbiased.pcorr_bound_up(isamp);
%         m_conf = pcorr_bound_up_given_choice(isamp);
%         m_choice_unbiased = 1;
%     case 2 %blo
%         m_like = probMapUnbiased.pcorr_bound_lo(isamp);
%         m_conf = pcorr_bound_lo_given_choice(isamp);
%         m_choice_unbiased = 0;
%     case 3 %notabs
%         m_like = probMapUnbiased.pcorr_nobound(idv,isamp);
%         % m_conf = pcorr_given_choice(idv,isamp);
%         m_conf = pcorr_given_choice(idv);
%         m_choice_unbiased = probMapUnbiased.best_choice_nobound(idv,isamp);
% end


%%


function m_like_right_separatrix = find_like_separatrix_given_belief_and_prior(K,prior_K,limit_prightbias)
% finds the separatrix between right and left bias reports 

m_like_right_separatrix = fminbnd(@eval_fun,0,1);


    function err = eval_fun(m_like_right)
        
        [~,pRightBias] = one_trial_bayesian_update(K,prior_K, 1 ,m_like_right);
        
        err = sqrt((pRightBias-limit_prightbias).^2);
        
    end
end

