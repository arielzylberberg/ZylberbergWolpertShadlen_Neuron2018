function h = errorbarxy(x,y,lerrx, uerrx, lerry, uerry)

n = length(x);
h = nan(n,2);
color = 'k';
for i=1:n
    hold all
    h(i,1) = plot([x(i)-lerrx(i),x(i)+uerrx(i)],[y(i),y(i)],'color',color);
    h(i,2) = plot([x(i),x(i)],[y(i)-lerry(i),y(i)+uerry(i)],'color',color);
end
