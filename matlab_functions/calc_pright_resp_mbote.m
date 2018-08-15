function right_conf_nobias = calc_pright_resp_mbote(dv,tt,drifts,pdrift)

if nargin<4 || isempty(pdrift)
    pdrift = ones(size(drifts));
end
pdrift = pdrift/sum(pdrift);

ppos = zeros(size(dv));
pall = zeros(size(dv));
for i=1:length(drifts)
    
    if drifts(i)>0
        ppos = ppos + pdrift(i) * normpdf(dv,drifts(i)*tt,sqrt(tt));
    elseif drifts(i)==0
        ppos = ppos + pdrift(i) * 1/2 * normpdf(dv,drifts(i)*tt,sqrt(tt));
    end
    
    pall = pall + pdrift(i) * normpdf(dv,drifts(i)*tt,sqrt(tt));
    
end

right_conf_nobias = ppos./pall;
