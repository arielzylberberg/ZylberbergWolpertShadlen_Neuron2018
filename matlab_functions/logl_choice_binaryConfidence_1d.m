function [logl,p_right_high,p_right_low,p_left_high,p_left_low] = ...
    logl_choice_binaryConfidence_1d(P,probMap,phi,choice,confidence_binary,dotdur,coh)
% function [logl,p_right_high,p_right_low,l_left_high,p_left_low] = logl_choice_binaryConfidence_1d(P,probMap,phi,choice,confidence_binary,dotdur,coh)

%% p choices - no sbet
% p_up = nan(size(choice));
% p_lo = nan(size(choice));
% p_out = nan(size(choice));
% I = sbet_on==0;
% [~,p_up(I),p_lo(I)] = logl_choice_VD_1d(P,choice(I),dotdur(I),coh(I));

%% p choices & conf
p  = P.notabs.pdf;
Up = P.up.pdf_t;
Lo = P.lo.pdf_t;
drift = P.drift;
y = P.y;

lastStep = ceil(dotdur/(P.t(2)-P.t(1)));

% I = sbet_on==1;
% [p_up(I),p_lo(I),p_out(I),p_up_t,p_lo_t,p_out_t,belief_c_t,belief_nc_t] = ...
%         pResp_withSbet_levelcurve(p,Up,Lo,y,drift,phi,[],coh(I),lastStep(I),probMap);

[p_right_high,p_right_low,p_left_high,p_left_low] = ...
                pResp_withBinaryConfidence(p,Up,Lo,y,drift,phi,[],coh,lastStep,probMap);
            
%% calc likelihood
p_right_high(p_right_high<eps)   = eps;
p_right_low(p_right_low<eps)     = eps;
p_left_high(p_left_high<eps)   = eps;
p_left_low(p_left_low<eps)     = eps;
conf = confidence_binary;

pPred = p_right_high.*(choice==1 & conf==1) + ...
        p_right_low.*(choice==1 & conf==0) + ...
        p_left_high.*(choice==0 & conf==1) + ...
        p_left_low.*(choice==0 & conf==0);


logl = -sum(log(pPred)); %logl

