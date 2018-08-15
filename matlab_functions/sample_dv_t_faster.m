function [idv,samp_num,CODE] =  sample_dv_t_faster(P,dotdur,icoh)
% function [idv,samp_num,CODE] =  sample_dv_t_faster(P,dotdur,icoh)
    
drifts = P.drift(icoh);
ntr = length(icoh);
nt = length(P.t);
dt = P.t(2)-P.t(1);
samp_num = ceil(dotdur/dt);
ev = bsxfun(@plus,randn(ntr,nt)*sqrt(dt),drifts*dt);
cev = cumsum(ev,2);

dy = P.y(2)-P.y(1);

Bup = P.Bup(1); % assumes flat bounds !!!
A = [cev>=(Bup-dy/2), true(ntr,1)];
[I,J] = find(A==1);
[~,m] = unique(I, 'first');
dt_upperbound = J(m);

Blo = P.Blo(1); % assumes flat bounds !!!
A = [cev<=(Blo+dy/2), true(ntr,1)];
[I,J] = find(A==1);
[~,m] = unique(I, 'first');
dt_lowerbound = J(m);


IND = sub2ind(size(cev),[1:ntr]',samp_num);
dv = cev(IND);
CODE = ones(ntr,1)*3;%notabs

I = dt_upperbound<dt_lowerbound & dt_upperbound<=samp_num;
CODE(I) = 1;
dv(I) = Bup;
samp_num(I) = dt_upperbound(I);

I = dt_upperbound>dt_lowerbound & dt_lowerbound<=samp_num;
CODE(I) = 2;
dv(I) = Blo;
samp_num(I) = dt_lowerbound(I);

idv = findclose(P.y,dv);

end