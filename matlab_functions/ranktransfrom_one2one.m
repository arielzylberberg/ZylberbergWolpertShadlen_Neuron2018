function yy = ranktransfrom_one2one(y,v)
% function yy = ranktransfrom_one2one(y,v)

yy = nan(size(y));
n = size(y,2);
[~,I] = sort(v);
for i=1:n
    [~,J] = sort(y(:,i));
    yy(J,i) = v(I);
end