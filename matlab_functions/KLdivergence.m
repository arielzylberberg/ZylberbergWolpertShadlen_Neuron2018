function KLdiv = KLdivergence(P,Q)
% P y Q son counts, no probabilidades !!

% sanity check
if all(P<1) || all(Q<1)
    error('los inputs deben ser counts');
end

% corrigo para que no haya ceros, agregando 1 al count
if any(P==0) || any(Q==0)
    P = P+1;
    P = P / sum(P); % renorm
    Q = Q+1;
    Q = Q / sum(Q); % renorm
end

% normalizo los counts a probabilidades
P = P / sum(P);
Q = Q / sum(Q); 

% calc
KLdiv = sum(log(P./Q).*P);