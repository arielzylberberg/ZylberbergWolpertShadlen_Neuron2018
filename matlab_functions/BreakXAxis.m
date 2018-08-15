function BreakXAxis(vhandle,center,varargin)
% function BreakXAxis(handle,center)
% handle: axis handle
% center: where to center the X-break, in axis coordenates

% if nargin<1
%     handle=gca;
%     center=0.021;
% elseif nargin<2
%     center = 0.021;
% end

break_scaling = 1;
for i=1:length(varargin)
    if isequal(varargin{i},'break_scaling')
        break_scaling = varargin{i+1};
    end
end

for i=1:length(vhandle)
    
    handle = vhandle(i);
    axes(handle);
    
    posit = get(handle,'position');
    w = posit(3);
    spacing = w/50; % in fig coords
    height = w/50 * break_scaling;
    yli = ylim;
    [xx,Yf] = ax2fig(handle,center,yli(1));
    
    
    Xf = calc_equal_spacing_log(xx,spacing);
    
    
    Y1 = [Yf - height/2;
        Yf + height/2];
    
    % Make the bar and ticks
    h1 = annotation(gcf, 'rectangle', [Xf(1), Y1(1), diff(Xf), diff(Y1)]);
    annotation(gcf, 'line', [Xf(1)-spacing/2 Xf(1)+spacing/2], [Y1(1) Y1(2)]);
    annotation(gcf, 'line', [Xf(2)-spacing/2 Xf(2)+spacing/2], [Y1(1) Y1(2)]);
    
    set(gcf,'color','w')
    set(h1,'facecolor','w','edgecolor','none')
end