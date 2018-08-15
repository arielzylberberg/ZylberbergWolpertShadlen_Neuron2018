function [conf_unbiased_right,conf_unbiased] = calc_like_from_conf(prior_K,K,choice,conf)

%calc like from conf

conf_right = choice.*conf+(1-choice).*(1-conf);

% for i=1:length(choice)
%     
%     pnextright = prior_K(i,:)*K(:);
%     ev_ratio = conf_right(i)/(1-conf_right(i))*(1-pnextright)/pnextright;
%     if ev_ratio>10000
%         conf_unbiased_right = 1;
%     else
%         conf_unbiased_right = ev_ratio/(ev_ratio+1);
%     end
% end

conf_unbiased_right = nan(size(choice));
for i=1:length(choice)
    
    pnextright = prior_K(i,:)*K(:);
    ev_ratio = conf_right(i)/(1-conf_right(i))*(1-pnextright)/pnextright;
    if ev_ratio>10000
        conf_unbiased_right(i) = 1;
    else
        conf_unbiased_right(i) = ev_ratio/(ev_ratio+1);
    end
end


conf_unbiased = choice.*conf_unbiased_right + (1-choice).*(1-conf_unbiased_right);


end