function bdata_modif = nonlin_mapping_belief(bdata,bmodel)
% bdata_modif = nonlin_mapping_belief(dat(1).belief,M(1).belief)
% inputs in the 0-1 range

fun = @ (theta) fun_error(theta,bdata,bmodel);

% bdata_modif = fn( 5 );
theta0 = 4;
options = optimset('Display','iter','FunValCheck','on');
theta = fminsearch(fun,theta0,options);

[~,bdata_modif] = fun(theta);

end

function [err,yp] = fun_error(theta,bdata,bmodel)

alpha = theta(1);
x = bdata;
y = bmodel;
yp = 1./(1+exp(-alpha*(x-0.5)));

dife = bsxfun(@minus,y,yp);
err = dife.^2;
err = nanmean(err(:));

end