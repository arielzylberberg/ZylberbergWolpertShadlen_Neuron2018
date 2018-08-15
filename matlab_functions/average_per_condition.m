function [Mean,Stdev,uni_conditions,tr_per_cond,idx_cond] = average_per_condition(data, conditions,varargin)
% function [Mean,Stdev,uni_conditions,tr_per_cond,idx_cond] = average_per_condition(data, conditions,varargin)
% dado un vector o matriz y una matriz de condiciones, devuelve
%media, desvio, etc. para cada condicion
% data: trials X time
% conditions: trials X conditionType
% varargin: {'filter',filtro}

tr_per_cond = []; % To-Do


filter = ones(size(data,1),1);
for i=1:length(varargin)
    if isequal(varargin{i},'filter')
        filter = varargin{i+1};
    end
end

globinds = filter==1;

% [uni_conditions,~,idx_cond_t] = unique(conditions(globinds,:),'rows');%before May 2017 !!
% the problem with this was the order of uni_conditions

% changed on May2017
if all(ismember(conditions(:),[0,1])) && all(sum(conditions,2)==1)
    n = size(conditions,2);
    idx_cond_t = nan(sum(globinds),1);
    for i=1:n
        I = conditions(globinds,i)==1;
        idx_cond_t(I) = i;
        uni_conditions = nan;
    end
else
    [uni_conditions,~,idx_cond_t] = unique(conditions(globinds,:),'rows');
    n = size(uni_conditions,1);
end


ntr = size(data,1);
nsamples = size(data,2);
% n = size(uni_conditions,1);

Mean  = nan(n,nsamples);
Stdev = nan(n,nsamples);

idx_cond = nan(ntr,1);
idx_cond(globinds) = idx_cond_t;
for i=1:n
   Mean(i,:)  = nanmean(data(idx_cond==i,:),1);
   Stdev(i,:) = nanstd(data(idx_cond==i,:),[],1);
end



