function [x_fig,y_fig]= ax2fig(handles,x,y)
% function [x_fig,y_fig]= ax2fig(handles,x,y)
% converts axes position to figure position
% Its useful for annotating graphs

pos=get(handles,'position');
xl=get(handles,'xlim');
yl=get(handles,'ylim');

x_fig=(x-xl(1))/(xl(2)-xl(1))*pos(3)+pos(1);
y_fig=(y-yl(1))/(yl(2)-yl(1))*pos(4)+pos(2);