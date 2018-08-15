function set_figure_size(hfig,xSize,filename)
% function set_figure_size(hfig,xSize,filename_optional)
% sets figure size, for printing. xSize in centimeters

dosave = true;
if nargin<3 || isempty(filename)
    dosave = false;
end

% MEDIDAS LETTER: xSize = 21.6; ySize = 27.9;
pos = get(hfig,'Position');
YoverX = pos(4)/pos(3);

set(hfig,'PaperUnits','centimeters','PaperPosition',[1 1 xSize YoverX*xSize])

if (dosave)
    saveas(hfig,filename,'epsc');
end
