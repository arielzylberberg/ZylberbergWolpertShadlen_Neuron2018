function li = same_lim(hax)

if nargin==0
    hax=gca;
end

xli = get(hax,'xlim');
yli = get(hax,'ylim');

li(1)=min([xli(:);yli(:)]);
li(2)=max([xli(:);yli(:)]);

xlim(li);
ylim(li);