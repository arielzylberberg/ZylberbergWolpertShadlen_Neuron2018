function probMap = beliefmap_1d_sdrift(P,Up,Lo,y,drift,prior)
% function probMap = beliefmap_1d_sdrift(P,Up,Lo,y,drift,prior)
% P [ncoh,ny,nt]
% Up [ncoh,nt]
% Lo [ncoh,nt]
% prior: prob. of each unique signed coherence
% 

%prior over coherences
prior = prior/nansum(prior);
signed = drift;

%at bound, up
aux = bsxfun(@times,Up,prior(:)); %p(<x,t>|coh).p(coh)
pcoh_xt_bound_up = bsxfun(@times,aux,1./sum(aux,1)); % p(coh|<x,t>), given x at bound

%same, lo
aux = bsxfun(@times,Lo,prior(:)); %p(<x,t>|coh).p(coh)
pcoh_xt_bound_lo = bsxfun(@times,aux,1./sum(aux,1));% p(coh|<x,t>), given x at bound

%probability of correct with bound reached at <x,t>
pcorr_bound_lo = sum(pcoh_xt_bound_lo(signed<0,:));
pcorr_bound_up = sum(pcoh_xt_bound_up(signed>0,:));
if any(signed==0)
    pcorr_bound_lo = pcorr_bound_lo + 0.5*sum(pcoh_xt_bound_lo(signed==0,:),1);
    pcorr_bound_up = pcorr_bound_up + 0.5*sum(pcoh_xt_bound_up(signed==0,:),1);
end

%when not at bound
aux = bsxfun(@times,P, reshape(prior,length(prior),1,1));
pcoh_xt_notbound = bsxfun(@times,aux,1./sum(aux,1));

[ntr,ny,nt] = size(P);
pcorr_nobound = nan(ny,nt);
pcorr_nobound(y>0,:) = squeeze(sum(pcoh_xt_notbound(signed>0,y>0,:)));
pcorr_nobound(y<=0,:) = squeeze(sum(pcoh_xt_notbound(signed<0,y<=0,:)));
if any(signed==0)
    pcorr_nobound  = pcorr_nobound + 0.5*squeeze(sum(pcoh_xt_notbound(signed==0,:,:),1));
end

pcorr_nobound = clip(pcorr_nobound,eps,1-eps);

% For unequal priors, the best 
% response might not be given by the sign of y.
% In this case, pcorr will be below 0.5. 
% Identify these cases, and invert the best response
best_choice_nobound = zeros(ny,nt);
%what we assumed so far
best_choice_nobound(y>0,:) = 1;
best_choice_nobound(y<=0,:) = 0;
inds = pcorr_nobound<0.5;
pcorr_nobound(inds) = 1 - pcorr_nobound(inds);%flip conf
best_choice_nobound(inds) = 1 - best_choice_nobound(inds);%flip choice


probMap = struct('pcorr_bound_lo',pcorr_bound_lo,'pcorr_bound_up',pcorr_bound_up,...
    'pcorr_nobound',pcorr_nobound,'best_choice_nobound',best_choice_nobound);

