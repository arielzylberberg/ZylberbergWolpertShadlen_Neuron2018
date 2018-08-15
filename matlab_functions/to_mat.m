function mat = to_mat(var,session,nses)
% function mat = to_mat(var,session,nses)

ntr = max(Rtable(session(~isnan(session))));
u = nanunique(session);
if nargin<3
    nses = length(u);
end

mat = nan(ntr,nses);
for i=1:length(u)
    inds = session==u(i);
    mat(1:sum(inds),i) = var(inds);
end



