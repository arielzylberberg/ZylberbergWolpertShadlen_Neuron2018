function [err,P,probMap,p_right_high,p_right_low,p_left_high,p_left_low,for_plotting] = ...
    wrapper_dtb_parametricbound_vd_choice_confidenceBinary(theta,dotdur,coh,choice,c,confidence_binary,prior,pars,plot_flag)

% dt_flag = pars.dt_flag;
% error_flag = pars.error_flag;

%%
kappa  = theta(1);
B0     = theta(2);
a      = theta(3);
d      = theta(4);
coh0   = theta(5);
y0a    = theta(6);
phi    = theta(7); % confidence separatrix in prob. scale

%%
notabs_flag = 1;

%%
dt = 0.0005;
t  = 0:dt:[nanmax(dotdur)+dt];
y  = linspace(-3,3,1500)';%esto puede ser problematico
y0 = zeros(size(y));
y0(findclose(y,0))=1;
y0 = y0/sum(y0);

%% bounds
[Bup,Blo] = expand_bounds(t,B0,a,d,'Exponential');

%%
% prior = Rtable(coh)/sum(Rtable(coh));

%%
% drift = kappa * unique(coh - coh0);
drift = kappa * unique(coh + coh0);
% notabs_flag = false;

% P = feval(F,drift,t,prior,hazard,y,y0,notabs_flag);
% P = dtb_fp_cn_searchbnd2(drift,t,prior,q,y,y0,notabs_flag);
% P = dtb_fp_cn_vec(drift,t,Bup,Blo,y,y0,notabs_flag);
P = dtb_fp_cc_vec(drift,t,Bup,Blo,y,y0,notabs_flag);
% P = spectral_dtbAA(drift,t,Bup,Blo,y,y0,notabs_flag);

%%
p = P.notabs.pdf;
Up = P.up.pdf_t;
Lo = P.lo.pdf_t;
drift = P.drift;
y = P.y;

probMap = beliefmap_1d_sdrift(p,Up,Lo,y,drift,prior); % get map
% [err,p_up,p_lo,p_out] = logl_choice_sbet_VD_1d(P,probMap,phi,choice,dotdur,coh,sbet_on,optedout);

[err,p_right_high,p_right_low,p_left_high,p_left_low] = ...
    logl_choice_binaryConfidence_1d(P,probMap,phi,choice,confidence_binary,dotdur,coh);

%% print
fprintf('err=%.3f kappa=%.2f B0=%.2f a=%.2f d=%.2f coh0=%.2f y0=%.2f phi=%.2f \n',...
    err,kappa,B0,a,d,coh0,y0a,phi);

%%
for_plotting = [];
if plot_flag
    
    [coh_fine,prh_fine,prl_fine,plh_fine,pll_fine] = evalmoredrifts_vd_confBinary(theta,dotdur,probMap);
    
    %%
    
    figure(1);clf
    
    subplot(1,3,1);
    color = get(0,'DefaultAxesColorOrder');
    curva_media(choice==1&confidence_binary==1,coh,[],3);
    hold all
    curva_media(choice==1&confidence_binary==0,coh,[],3);
    curva_media(choice==0&confidence_binary==1,coh,[],3);
    curva_media(choice==0&confidence_binary==0,coh,[],3);
    hold all
    plot(coh_fine,prh_fine,'-','color',color(1,:))
    plot(coh_fine,prl_fine,'-','color',color(2,:))
    plot(coh_fine,plh_fine,'-','color',color(3,:))
    plot(coh_fine,pll_fine,'-','color',color(4,:))
    xlim([-0.55,0.55])
    
    %
    %     coh_side = coh_fine(coh_fine>=0);
    %     pc_fine = nan(size(coh_fine));
    %     pc_fine(coh_fine>0) = prh_fine(coh_fine>0) + prl_fine(coh_fine>0);
    %     pc_fine(coh_fine<0) = plh_fine(coh_fine<0) + pll_fine(coh_fine<0);
    %     [~,pc_fine] = curva_media(pc_fine,abs(coh_fine),[],1);
    
    subplot(1,3,2);
    p_correct_rightward = double(coh_fine>0);
    p_correct_rightward(coh_fine==0) = 0.5;
    [coh_unsigned, pc_unsigned, phigh_given_c_unsigned, phigh_given_nc_unsigned] = ...
        split_into_acc_conf(coh_fine,prh_fine,prl_fine,plh_fine,pll_fine,p_correct_rightward);
    
    
    for_plotting = struct('coh_unsigned',coh_unsigned, 'pc_unsigned', pc_unsigned,...
        'phigh_given_c_unsigned',phigh_given_c_unsigned, 'phigh_given_nc_unsigned',phigh_given_nc_unsigned);
    
    correct = (choice==1).*(coh>0)+(choice==0).*(coh<0);
    correct(coh==0) = 0.5;
    
%     curva_media(correct,abs(coh),[],3);
%     hold all
%     plot(coh_unsigned,pc_unsigned,'b');
%     yy = pc_unsigned(coh_unsigned==0);
%     plot([0.016,0.019],[yy,yy],'b');
    
    dotsanalysis.plot_log(coh,correct,coh_unsigned,pc_unsigned,'break_scaling',4);
    title('accuracy')
    ylabel('p correct')
    
    subplot(1,3,3);
    
    I = correct==0 | coh==0;
    dotsanalysis.plot_log(coh(I),confidence_binary(I),coh_unsigned,phigh_given_nc_unsigned,'break_scaling',4,'color','r')
    I = correct==1 | coh==0;
    dotsanalysis.plot_log(coh(I),confidence_binary(I),coh_unsigned,phigh_given_c_unsigned,'break_scaling',4,'color','b')

%     coh_plot=coh;
%     coh_plot(coh==0)=0.016;%for log scale
%     curva_media(confidence_binary==1,abs(coh_plot),correct==1 | coh==0,3);
%     hold all
%     plot(coh_unsigned,phigh_given_c_unsigned,'b');
%     hold on
%     yy = phigh_given_c_unsigned(coh_unsigned==0);
%     plot([0.016,0.019],[yy,yy],'b');
%     
%     curva_media(confidence_binary==1,abs(coh_plot),correct==0 & coh~=0,3);
%     hold all
%     plot(coh_unsigned,phigh_given_nc_unsigned,'r');
    title('confidence')
    ylabel('p high')
    
    set(gcf,'Position',[270   775  1135   248])
    format_figure(gcf,'FontSize',15);
    
    ch = get(gcf,'children');
    set(ch(1:2),'xscale','log');
    drawnow
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coh_unsigned, pc_unsigned, phigh_given_c_unsigned, phigh_given_nc_unsigned] = ...
    split_into_acc_conf(coh_signed,prh,prl,plh,pll,p_correct_rightward)
% function [coh_side, pc_unsigned, phigh_given_c_unsigned, phigh_given_nc_unsigned] = ...
%     split_into_acc_conf(coh_signed,prh,prl,plh,pll,p_correct_rightward)
% converts the proportion over choice & confidence into accuracy,
% phigh_correct and phigh_incorrect


coh_signed = coh_signed(:);
p_correct_rightward = p_correct_rightward(:);


coh_unsigned = coh_signed(coh_signed>=0);
pcr = p_correct_rightward; % for each element in coh_signed, prob that the req resp is rightward


pc = prh.*pcr + prl.*pcr + plh.*(1-pcr) + pll.*(1-pcr);
[~,pc_unsigned] = curva_media(pc,abs(coh_signed),[], 0 );


phigh_given_c = pcr.*prh./(prh+prl)+(1-pcr).*plh./(plh+pll);
phigh_given_nc = (1-pcr).*prh./(prh+prl)+(pcr).*plh./(plh+pll);
% phigh_given_c = (prh.*pcr+plh.*(1-pcr))./(prh.*pcr+plh.*(1-pcr)+prl.*pcr+pll.*(1-pcr));
% phigh_given_nc = (prh.*(1-pcr)+plh.*pcr)./(prh.*pcr+plh.*(1-pcr)+prl.*pcr+pll.*(1-pcr));

[~,phigh_given_c_unsigned] = curva_media(phigh_given_c,abs(coh_signed),[], 0 );
[~,phigh_given_nc_unsigned] = curva_media(phigh_given_nc,abs(coh_signed),[], 0 );

end