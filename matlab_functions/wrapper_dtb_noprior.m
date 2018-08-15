function P = wrapper_dtb_noprior(theta,dotdur,coh,signedzero_flag)

%%
kappa  = theta(1);
B0     = theta(2);
coh0   = theta(3);
y0a    = theta(4);

%%
dt = 0.0005;
maxtime = nanmax(dotdur)+0.1;
notabs_flag = 1;
t  = 0:dt:[maxtime + dt];
Bup = B0*ones(size(t));
Blo = -1*Bup;
y  = symmetric_scale(max(Bup)*1.2,0.004);
y0 = zeros(size(y));
y0(findclose(y, y0a )) = 1;
y0 = y0/sum(y0);

if signedzero_flag
    uacoh = nanunique(abs(coh));
    uacoh(uacoh==0) = 10^-9;
    ucoh = sort([-uacoh(:);uacoh(:)]);
else
    ucoh = nanunique(coh);
end
drift = kappa * unique(ucoh + coh0);
P = dtb_fp_cc_vec(drift,t,Bup,Blo,y,y0,notabs_flag); % could be replaced by analytic dtb

