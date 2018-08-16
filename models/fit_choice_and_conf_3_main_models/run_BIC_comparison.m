clear

%%
modeldir = './';
% folders = {'use_conf_as_like_nonparW',...
%     'ignore_conf_nonparW',...
%     'use_conf_as_like_flatW',...
%     'use_conf_as_like_gaussW',...
%     'use_conf_as_like_singleW',...
%     'use_conf_as_like_twoW',...
%     'sample_one_like_nonparW'};

% folders = {'ignore_conf_nonparW_shrink',...
%            'use_like_nonparW_shrink',...
%            'use_conf_nonparW_shrink',...
%            'use_like_nonparW_shrink_extremes',...
%            'use_rand_nonparW_shrink',...
% };

folders = {'ignore_conf_nonparW_shrink',...
           'use_like_nonparW_shrink',...
           'use_conf_nonparW_shrink',...
};

numParam = [5,5,5];

include_models = [1:length(numParam)];
% include_models = [2,4,6,7,8];


folders = folders(include_models);
numParam = numParam(include_models);

aux = load('../../data/data_CCONF.mat');
% aux = load('../../models5_clean/prep_data_for_model/data_for_model.mat');
numObs = Rtable(cat(1,aux.dat.group));
    
nsujs = length(numObs);
for i=1:length(folders)
    for j=1:nsujs
        dat = load(fullfile(modeldir,folders{i},['suj_',num2str(j)]),'fval','theta');
        logl(i,j) = -1 * dat.fval;
        theta{i,j} = dat.theta;
    end
end


Ibase = 1;

for j=1:nsujs
    [aic_base(j),bic_base(j)] = aicbic(logl(Ibase,j),numParam(Ibase),numObs(j));
    for i=1:length(folders)
        [aic(i,j),bic(i,j)] = aicbic(logl(i,j),numParam(i),numObs(j));
    end
end

disp('diff AIC:')
aic_diff = bsxfun(@plus,aic,-aic_base)
disp('diff BIC:')
bic_diff = bsxfun(@plus,bic,-bic_base)

disp('log')
best = bsxfun(@eq,logl,max(logl))
v = -1 * round(bsxfun(@plus,logl,-logl(best)'));
T = table(v(:,1),v(:,2),v(:,3),'RowNames',folders,'VariableNames',{'S1','S2','S3'})

disp('best per suj:')
best = bsxfun(@eq,bic,min(bic))

v = round(bsxfun(@plus,bic,-bic(best)'));
T = table(v(:,1),v(:,2),v(:,3),'RowNames',folders,'VariableNames',{'S1','S2','S3'})


%% 

params_bayesian_model = cat(1,theta{2,:})'

