function [p,T,X,S] = curva_media_wrapper(DEP,INDEP,filt_mat,plotflag)
%function [p,T,X,S] = curva_media_wrapper(DEP,INDEP,filt_mat,plotflag)


if isvector(filt_mat)
    filt_mat = adummyvar(filt_mat,1)==1;
end

for i=1:size(filt_mat,2)
    inds = filt_mat(:,i);
    [T{i},X{i},S{i}] = curva_media(DEP,INDEP,inds,0);
end

if plotflag
    p = publish_plot(1,1);
else
    p = [];
end
if plotflag
    colores = rainbow_colors(length(T),'colorType',0);
    for i=1:length(T)
        if plotflag==1
            plot(T{i},X{i},'.-','color',colores(i,:))
        elseif plotflag==2
            errorbar(T{i},X{i},S{i},'.-','color',colores(i,:))
        end
        hold on
    end
end

if plotflag
    p.format('FontSize',15);
end