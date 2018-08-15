function right_like = calc_pright_resp(idv,samp_num,CODE,probMap)
% function right_like = calc_pright_resp(idv,samp_num,CODE,probMap)
% likehood of a rightward response given the dv and time, using the map

[ntr,nsamp] = size(idv);
idv      = idv(:);
samp_num = samp_num(:);
CODE = CODE(:);

% for the ones which toched the bound, used the idv that is just
% below/above depending on the bound
idv(CODE==1) = idv(CODE==1) - 1; % upper bound
idv(CODE==2) = idv(CODE==2) + 1; % lower bound

[fil,col] = size(probMap.pcorr_given_right_nobound);

IND = sub2ind([fil,col],idv,samp_num);

right_like = probMap.pcorr_given_right_nobound(IND);

%resize
right_like = reshape(right_like,[ntr,nsamp]);

end