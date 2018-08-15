function [deltab,belief_shifted] = calc_deltab(belief,session)
% [deltab,belief_shifted] = calc_deltab(belief,session)

base = 0.5;
% if any(belief<0)
%      base = 0;
% end

deltab = nan(size(belief));
belief_shifted = nan(size(belief));

I = ~isnan(belief);
s = session(I);
b = belief(I);
uni = nanunique(s);
nses = length(nanunique(s));
db = nan(size(b));
b_shifted = nan(size(b));
for i=1:nses
    inds = find(s==uni(i));
    db(inds) = diff([base; b(inds)]);
    
    b_shifted(inds) = [base; b(inds(1:end-1))];
    
end
deltab(I) = db;
belief_shifted(I) = b_shifted;