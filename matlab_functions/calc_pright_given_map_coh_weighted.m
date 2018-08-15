function p_right = calc_pright_given_map_coh_weighted(probMap,P,icoh,prior_K_prev,it)
% function p_right_t = calc_pright_given_map_coh(probMap,P,icoh,imap)
% computes the prob of a rightward choice given the coherence index, 
% marginalizing over the evidence and world states.


pdf = squeeze(P.notabs.pdf(icoh,:,it))';

n = length(probMap);
prior_K_prev = prior_K_prev/sum(prior_K_prev);

pright_nobound = zeros(size(pdf));
pright_bound_up = zeros(1,it);
pright_bound_lo = zeros(1,it);

pleft_nobound = zeros(size(pdf));
pleft_bound_up = zeros(1,it);
pleft_bound_lo = zeros(1,it);

for i=1:n
    pright_nobound = pright_nobound + prior_K_prev(i) * probMap(i).pev_notabs(:,it) .* probMap(i).pright_nobound(:,it);
    pright_bound_up = pright_bound_up + prior_K_prev(i) * probMap(i).pev_bound_up(1:it) .* probMap(i).pright_bound_up(1:it);
    pright_bound_lo = pright_bound_lo + prior_K_prev(i) * probMap(i).pev_bound_lo(1:it) .* probMap(i).pright_bound_lo(1:it);
    
    pleft_nobound = pleft_nobound + prior_K_prev(i) * probMap(i).pev_notabs(:,it) .* (1-probMap(i).pright_nobound(:,it));
    pleft_bound_up = pleft_bound_up + prior_K_prev(i) * probMap(i).pev_bound_up(1:it) .* (1-probMap(i).pright_bound_up(1:it));
    pleft_bound_lo = pleft_bound_lo + prior_K_prev(i) * probMap(i).pev_bound_lo(1:it) .* (1-probMap(i).pright_bound_lo(1:it));
end

pright_nobound = pright_nobound./(pright_nobound+pleft_nobound);
pright_bound_up = pright_bound_up./(pright_bound_up+pleft_bound_up);
pright_bound_lo = pright_bound_lo./(pright_bound_lo+pleft_bound_lo);

    
ch_nobound = pright_nobound>0.5;

% not absorbed
right_notabs = sum(ch_nobound.*pdf);

% absorbed at bound
right_abs = sum ( P.up.pdf_t(icoh,1:it).*(pright_bound_up>0.5) + ...
    P.lo.pdf_t(icoh,1:it).*(pright_bound_lo>0.5) );

% sum
p_right = right_notabs + right_abs;