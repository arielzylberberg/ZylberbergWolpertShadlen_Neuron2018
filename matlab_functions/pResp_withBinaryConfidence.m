function [p_right_high,p_right_low,p_left_high,p_left_low] = ...
                pResp_withBinaryConfidence(P,Up,Lo,y,drift,phi,prior,coh,lastStep,probMap)
            
% function pResp_withSbet_levelcurve(P,Up,Lo,y,drift,phi,p_uscoh,coh,lastStep,probMap)
% P [ncoh,ny,nt]
% Up [ncoh,nt]
% Lo [ncoh,nt]
% p_uscoh prob. of each unique signed coherence
% lastStep: tdots/samp_dur
%

if nargin<10
    probMap = [];
end

[nd,ny,nt] = size(P);

if isempty(probMap)
    probMap = beliefmap_1d_sdrift(p,Up,Lo,y,drift,prior);
    
end
pcorr_bound_lo = probMap.pcorr_bound_lo(1:nt);
pcorr_bound_up = probMap.pcorr_bound_up(1:nt);
pcorr_nobound  = probMap.pcorr_nobound(:,1:nt);
best_choice_nobound = probMap.best_choice_nobound(:,1:nt);



% collapse based on the value of phi
pRight_high = zeros(ny,nt);
pRight_high(pcorr_nobound>phi & best_choice_nobound==1) = 1;

pLeft_high = zeros(ny,nt);
pLeft_high(pcorr_nobound>phi & best_choice_nobound==0) = 1;

pRight_low = zeros(ny,nt);
pRight_low(pcorr_nobound<phi & best_choice_nobound==1) = 1;

pLeft_low = zeros(ny,nt);
pLeft_low(pcorr_nobound<phi & best_choice_nobound==0) = 1;


pRightBound_high = zeros(nt,1);
pLeftBound_high = zeros(nt,1);
pRightBound_low = zeros(nt,1);
pLeftBound_low = zeros(nt,1);

pRightBound_high(pcorr_bound_up>phi) = 1;
pLeftBound_high(pcorr_bound_lo>phi) = 1;
pRightBound_low(pcorr_bound_up<phi) = 1;
pLeftBound_low(pcorr_bound_lo<phi) = 1;

%now evaluate each drift and time
p_up_t_high = cumsum(bsxfun(@times,Up,pRightBound_high'),2) + squeeze(sum(bsxfun(@times,P,reshape(pRight_high,[1,size(pRight_high)])),2));
p_lo_t_high = cumsum(bsxfun(@times,Lo,pLeftBound_high'),2) + squeeze(sum(bsxfun(@times,P,reshape(pLeft_high,[1,size(pLeft_high)])),2));

p_up_t_low = cumsum(bsxfun(@times,Up,pRightBound_low'),2) + squeeze(sum(bsxfun(@times,P,reshape(pRight_low,[1,size(pRight_low)])),2));
p_lo_t_low = cumsum(bsxfun(@times,Lo,pLeftBound_low'),2) + squeeze(sum(bsxfun(@times,P,reshape(pLeft_low,[1,size(pLeft_low)])),2));

% p_out_t = squeeze(sum(bsxfun(@times,P,reshape((1.0-pUp-pLo),[1,size(pUp)])),2)); %how it was before May-2016
% p_out_t = cumsum(bsxfun(@times,Up,(1-pRightBound_high)'),2) + cumsum(bsxfun(@times,Lo,(1-pLeftBound_high)'),2) + ...
%     squeeze(sum(bsxfun(@times,P,reshape((1.0-pRight_high-pLeft_high),[1,size(pRight_high)])),2));

%
if isvector(drift)
    [ucoh,~,icoh] = unique(coh);
else
    icoh = 1:nd;
end

% % subj. prob. correct conditional on
% % coh and choice
% cUp = cumsum(Up,2);
% Iup = pcorr_nobound>phi & probMap.best_choice_nobound==1;
% norm_up = cUp + squeeze(sum(bsxfun(@times,P,reshape(Iup,[1,ny,nt])),2));
% aux = bsxfun(@times,Up,pcorr_bound_up);
% aux(isnan(aux))=0;
% pcbup = cumsum(aux,2);
% aux = reshape(pcorr_nobound.*Iup,[1,ny,nt]);
% aux(isnan(aux)) = 0;
% pc_up_t = pcbup + squeeze(sum(bsxfun(@times,P,aux),2));
% pc_up_t = pc_up_t./norm_up;
% 
% %same, lower
% cLo = cumsum(Lo,2);
% Ilo = pcorr_nobound>phi & probMap.best_choice_nobound==0;
% norm_lo = cLo + squeeze(sum(bsxfun(@times,P,reshape(Ilo,[1,ny,nt])),2));
% aux = bsxfun(@times,Lo,pcorr_bound_lo);
% aux(isnan(aux))=0;
% pcblo = cumsum(aux,2);
% aux = reshape(pcorr_nobound.*Ilo,[1,ny,nt]);
% aux(isnan(aux)) = 0;
% pc_lo_t = pcblo + squeeze(sum(bsxfun(@times,P,aux),2));
% pc_lo_t = pc_lo_t./norm_lo;
% 
% % transform up y lo to corr and ncorr
% belief_c_t = bsxfun(@times,pc_up_t,ucoh>0) + bsxfun(@times,pc_lo_t,ucoh<0);
% belief_nc_t = bsxfun(@times,pc_up_t,ucoh<0) + bsxfun(@times,pc_lo_t,ucoh>0);
%     
% 
% if sum(ucoh==0)>0
%     belief_c_t = belief_c_t + 0.5 * bsxfun(@times,pc_up_t(ucoh==0,:) + pc_lo_t(ucoh==0,:),ucoh==0);
%     belief_nc_t = belief_nc_t + 0.5 * bsxfun(@times,pc_up_t(ucoh==0,:) + pc_lo_t(ucoh==0,:),ucoh==0);
% end

%% now find for lastStep

ind = sub2ind(size(p_up_t_high),icoh,lastStep);
p_right_high = p_up_t_high(ind);
p_left_high = p_lo_t_high(ind);
p_right_low = p_up_t_low(ind);
p_left_low = p_lo_t_low(ind);
% p_out = p_out_t(ind);
% belief_c = belief_c_t(ind);
% belief_nc = belief_nc_t(ind);


