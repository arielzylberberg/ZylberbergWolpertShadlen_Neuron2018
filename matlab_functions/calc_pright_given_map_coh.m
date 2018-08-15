function p_right_t = calc_pright_given_map_coh(probMap,P,icoh)
% function p_right_t = calc_pright_given_map_coh(probMap,P,icoh,imap)
% computes the prob of a rightward choice given the coherence index
% and the confidence map and decision policy based on the map

% icoh = 4;
% imap = 41;

pdf = squeeze(P.notabs.pdf(icoh,:,:));
ch_nobound = probMap.best_choice_nobound;

% not absorbed
right_notabs = sum(ch_nobound.*pdf);

% absorbed at bound
right_abs = cumsum ( P.up.pdf_t(icoh,:).*(probMap.pcorr_bound_up>0.5) + ...
    P.lo.pdf_t(icoh,:).*(probMap.pcorr_bound_lo<0.5) );

% sum
p_right_t = right_notabs + right_abs;