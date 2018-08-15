function [h_ax,hx,hy] = nontouching_spines(h_axes,varargin)
% nontouching_spines(h_ax,varargin)

separation = 0.02;
extra = 0.0; %para que se vea el texto
do_X = true; % displace x axis
do_Y = true; % displace y axis
do_shrink = true; % shrink current axis
tl = 0.03;
for i=1:length(varargin)
    if isequal(varargin{i},'separation')
        separation = varargin{i+1};
    elseif isequal(varargin{i},'extra')
        extra = varargin{i+1};
    elseif isequal(varargin{i},'only_X')
        do_Y = false;
    elseif isequal(varargin{i},'only_Y')
        do_X = false;
    elseif isequal(varargin{i},'dont_shrink')    
        do_shrink = false;
    elseif isequal(lower(varargin{i}),'ticklength')    
        tl = varargin{i+1};
    end
end
% plot([-1.1 1],[0 1]);%,'LineSmoothing','on');


for i=1:length(h_axes)
    h_ax = h_axes(i);
    set(gcf,'CurrentAxes',h_ax);
    
    % fix axes lims
    set(gca,'XLimMode','manual')
    set(gca,'YLimMode','manual')
    
    pos = get(h_ax,'position');
    
    if do_shrink
        pos(1)=pos(1)+extra;
        pos(3)=pos(3)-extra;
        pos(2)=pos(2)+extra;
        pos(4)=pos(4)-extra;
    end
    set(gca,'position',pos);
end

for i=1:length(h_axes)
    h_ax = h_axes(i);
    set(gcf,'CurrentAxes',h_ax);
    XLI = xlim;
    YLI = ylim;
    %     set(h_ax,'XLimMode','manual')
    %     set(h_ax,'YLimMode','manual')
    
    XLAB = get(get(gca,'xlabel'),'String');
    YLAB = get(get(gca,'ylabel'),'String');
    
    posplot = get(h_ax,'position');
    if do_shrink
        posplot(1) = posplot(1)+separation;
        posplot(2) = posplot(2)+separation;
        posplot(3) = posplot(3)-separation;
        posplot(4) = posplot(4)-separation;
    end
    set(gca,'position',posplot);
    
    xcolor_original = get(h_ax,'xcolor');
    %     ycolor_original = get(h_ax,'ycolor');
    FontSize = get(h_ax,'FontSize');
    
    props = get(h_ax);
    
    
    if verLessThan('matlab','8.4')
        to_copy_x = {'FontSize','XColor','XTick','XTickLabel','XScale'};
        to_copy_y = {'FontSize','YColor','YTick','YTickLabel','YAxisLocation'};
    else
        to_copy_x = {'FontSize','XColor','XTick','XTickLabel','XTickLabelRotation','XScale'};
        to_copy_y = {'FontSize','YColor','YTick','YTickLabel','YAxisLocation','YTickLabelRotation'};

    end
    
    %% aux, calc tick length
    pos = get(gca,'position');
    posfig = get(gcf,'position');
    ratio_wh = pos(3)/pos(4) * posfig(3)/posfig(4) * 1.2;
    
    
    %% x
    if do_X
        hx = axes; 			% Create an axes object in the figure
        pos = posplot;
        pos(2) = pos(2)-separation;
        pos(4) = 0.01;
        set(gca,'position',pos)
        %     set(gca,'ycolor','w','ytick',[],'tickdir','out','ticklength',[0.02 0.02],'xlim',XLI)
        set(gca,'ycolor','w','ytick',[],'tickdir','out','ticklength',[tl tl],'xlim',XLI)
        set(gca,'xcolor',xcolor_original,'color','none','FontSize',FontSize);
        for ii=1:length(to_copy_x)
            f = to_copy_x{ii};
            set(gca,f,props.(f));
        end
        %     set(gca,'FontAngle','italic')
        xla = xlabel(XLAB);
        set(xla,'FontAngle','normal')
        
        set(h_ax,'xcolor','w','xtick',[])
    end
    
    
    
    %% y
    if do_Y
        hy = axes;
        pos = posplot;
        if isequal(get(h_ax,'YAxisLocation'),'right')
            pos(1) = pos(1) + pos(3)+separation;
        else
            pos(1) = pos(1) - separation;
        end
        pos(3) = 0.01;
        set(gca,'position',pos)
        %     set(gca,'xcolor','w','xtick',[],'tickdir','out','ticklength',[0.02 0.02],'ylim',YLI)
        set(gca,'xcolor','w','xtick',[],'tickdir','out','ticklength',[tl,tl]*ratio_wh,'ylim',YLI)
        %     set(gca,'ycolor',ycolor_original,'color','none','FontSize',FontSize);
        set(gca,'color','none');
        for ii=1:length(to_copy_y)
            f = to_copy_y{ii};
            set(gca,f,props.(f));
        end
        
        set(get(gca,'YLabel'),'Rotation', get(get(h_ax,'YLabel'),'Rotation') )
        
        %     set(gca,'FontAngle','italic')
        
        yla = ylabel(YLAB);
        set(yla,'FontAngle','normal');
        
        %         set(h_ax,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
        set(h_ax,'ycolor','w','ytick',[])
    end
    
    %%
    
    set(h_ax,'xlim',XLI,'ylim',YLI)
    
    %     set(h_ax,'visible','off')
    
    
    set(gcf,'CurrentAxes',h_axes(i));
    
    
    
end

end


