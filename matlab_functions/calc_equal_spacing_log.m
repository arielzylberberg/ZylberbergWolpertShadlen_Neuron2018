function X = calc_equal_spacing_log(x,delta)

xright = x + delta;
xleft = exp(2*log(x)-log(xright));

X = [xleft,xright];