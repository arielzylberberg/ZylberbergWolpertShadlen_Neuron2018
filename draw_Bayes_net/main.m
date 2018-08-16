b = bayes_net_draw();
b.add_node(2,3,'B');
b.add_node(2,2,'d');
b.add_node(1,2,'s');
b.add_node(1.5,1,'e',1);
b.add_connection('B','d');
b.add_connection('d','e');
b.add_connection('s','e');


b.add_hyperparam(2,3.7,'k_B');
b.add_connection('k_B','B');

b.add_hyperparam(.2,2.0,'k_c');
b.add_connection('k_c','s');


b.add_text('B','$B$',-1.5);
b.add_text('d','$d$',-1.5);
b.add_text('s','$|c|$',-1.5);
b.add_text('e','$x$',-1.5);

b.add_text('k_B','$k_B$',-1.5);
b.add_text('k_c','$k_c$');



b.add_plate([1,1,1,1.7],0.6,'BLOCKS');
b.add_plate([1,1,1,1],0.55,'TRIALS');
b.add_plate([1,1,1,.1],0.5,'SAMPLES');
% b.add_plate([1,1,.5,.5],0.4);
b.format();

% add text to the right
yy = linspace(3.1,.5,10);
text(3,yy(2),['$B\sim$Uniform$(k_B)$'],'interpreter','latex','horizontalalignment','left','FontSize',16)
text(3,yy(4),'$p(d=1|B) = B$','interpreter','latex','horizontalalignment','left','FontSize',16)
text(3,yy(6),['$|c|\sim$Uniform$(k_c)$'],'interpreter','latex','horizontalalignment','left','FontSize',16)
text(3,yy(8),'$T\sim$TruncExp$(k_T) $','interpreter','latex','horizontalalignment','left','FontSize',16)
text(3,yy(10),'$x\sim \mathcal{N}(\kappa c \Delta t, \sqrt{\Delta t})$','interpreter','latex','horizontalalignment','left','FontSize',16)


str = {'Block bias: ','Motion direction:','Motion strength:','Motion duration:','Momentary evidence:'};
for i=1:length(str)
    text(2.7,yy(2*i-1),str{i},'horizontalalignment','left','FontSize',14,'color','r','interpreter','tex')
end

%%
xlim([-1,10]);
set(gcf,'Position',[142    83  1130   722])
axis off
