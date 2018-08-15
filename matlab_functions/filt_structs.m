function [all_filt,pers_filt] = filt_structs(all,pers,filt)

I = pers.ref_to_all(:,filt);
I = I(~isnan(I));

f = fieldnames(all);
for i=1:length(f)
    all_filt.(f{i}) = all.(f{i})(I);
end

if ~isempty(pers)
%     pfilt = false(size(pers.ref_to_all));
%     inds = ~isnan(pers.ref_to_all);
%     pfilt(inds) = filt(pers.ref_to_all(inds));
%     pfilt = reshape(pfilt,size(pers.group));
% 
%     if any(mean(pfilt)==0 | mean(pfilt)==1)==0
%         error('filter is not correctly specified for pers')
%     end
%     pfilt = pfilt(1,:);%first row

    f = fieldnames(pers);
    for i=1:length(f)
        pers_filt.(f{i}) = pers.(f{i})(:,filt);
    end
else
    pers_filt = [];
end



