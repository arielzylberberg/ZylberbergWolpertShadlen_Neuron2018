function [beta,idx,stats,x,LRT] = f_regression(depvar,boolNormalize,indepvar,testSignificance)
% [beta,idx,stats,x,LRT] = f_regression(depvar,boolNormalize,indepvar,testSignificance)

% sanity check
if isempty(boolNormalize)
    boolNormalize = zeros(length(indepvar)/2,1);
end

if not(length(boolNormalize)==(length(indepvar)/2))
    error('input sizes dont match')
end

if nargin <4 || isempty(testSignificance)
    testSignificance = [];
end

% parse independent variables
X       = [];
ini     = 1;
cont    = 0;
for i=1:2:length(indepvar)
    cont = cont + 1;
    tempvar             = indepvar{i+1};
    X                   = [X, tempvar];
    names{cont}         = indepvar{i};
    nCol                = size(tempvar,2);
    idx.(indepvar{i})   = [ini:(ini+nCol-1)];
    ini                 = ini + nCol;
end


% preprocessing of inputs
media = nanmean(X);
stdev = nanstd(X,[],1);

if all(boolNormalize==0)
    indepvars = X;
else
    Xnorm = mresta(X,media);
    Xnorm = mmulti(Xnorm,stdev);
    
    for i=1:length(names)
        inds = idx.(names{i});
        if boolNormalize(i) == 1
            indepvars(:,inds) = Xnorm(:,inds);
        else
            indepvars(:,inds) = X(:,inds);
        end
    end
end

constant = 'off';
estdisp = 'on';
if length(nanunique(depvar))==2
    logistic_flag = 1;
    disp('logistic regression')
    distr = 'binomial';
    link = 'logit';
    
else
    logistic_flag = 0;
    disp('linear regression')
    distr = 'normal';
    link = 'identity';
end

[beta,dev,stats] = glmfit(indepvars,depvar,distr,'link',link,'constant',constant,'estdisp',estdisp);
stats.DEV = dev;
x = indepvars;


%% LIKELIHOOD RATIO TESTS - only for logistic
LRT = [];
if not(isempty(testSignificance)) && logistic_flag
    n = 1; % just to indicate binary y
    yfit  = glmval(beta,indepvars,link,'constant',constant);
    logLikelihoodFullModel = nansum(log(binopdf(depvar,n,yfit)));
    %         logLikelihoodFullModel = logl_linear_model(depvar,yfit);
    
    varsTest = testSignificance.vars;
    % restricted models
    for i=1:length(varsTest)
        id = varsTest(i);
        xrestricted = x;
        xrestricted(:,idx.(names{id})) = [];
        
        [betares] = glmfit(xrestricted,depvar,distr,'link',link,'constant',constant,'estdisp',estdisp);
        yfit = glmval(betares,xrestricted,'logit','constant','off');
        logLikelihoodRestrictedModel = nansum(log(binopdf(depvar,n,yfit)));
        %         if (logistic_flag)
        %             logLikelihoodRestrictedModel = nansum(log(binopdf(depvar,n,yfit)));
        %         else
        %             logLikelihoodRestrictedModel = logl_linear_model(depvar,yfit);
        %         end
        
        % Likelihood ratio test
        LRT(i).variableName    = names{id};
        LRT(i).cstat           = -2*logLikelihoodRestrictedModel-(-2)*logLikelihoodFullModel;
        LRT(i).dfe             = length(idx.(names{id}));
        LRT(i).p                = 1 - chi2cdf(LRT(i).cstat,LRT(i).dfe);
        LRT(i).beta             = beta(idx.(names{id}));
        
    end
end
end

%% %%%%%%%%%%%%%%%%%%
%%% HELPER FN %%%%%%
%%%%%%%%%%%%%%%%%%%%

% STILL NOT IMPLEMENTED
function L = logl_linear_model(y,yfit)
N = length(y);
Residuals = y-yfit;
sigma = nanstd(Residuals);    %find std.dev. of residuals
%Calculate L using simple expression
L = -N*.5*log(2*pi) - N*log(sigma) - (1/(2*sigma.^2))*sum(Residuals.^2);

%set K covariance matrix
%     K_matrix =  eye(N) * sigma^2;
%     L = -N*.5*log(2*pi) - sum(log(diag(chol(K_matrix)))) ...
%            - .5*Residuals' / (K_matrix)* Residuals;

end

