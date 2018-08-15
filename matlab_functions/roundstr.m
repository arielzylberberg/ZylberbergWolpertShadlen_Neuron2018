function str = roundstr(x,ndecimals)
%function str = roundstr(x,ndecimals)
% rounds to ndecimals and converts to string
% useful for plotting

xround = round(x*10^ndecimals)/10^ndecimals;
str = num2str(xround);