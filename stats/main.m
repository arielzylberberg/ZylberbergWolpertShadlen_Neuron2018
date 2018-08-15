
addpath(genpath('../matlab_functions/'));

%%
uni_sujs = 1:3;
data_path_to_file = '../data/data_CCONF.mat';
[d] = get_data(uni_sujs,data_path_to_file);
struct2vars(d);

%% Figure 2A
% Logistic regression

filt = true(size(choice)); 

ntr = sum(filt);
depvar = choice(filt);
dummy_prior = adummyvar(block_prior);
dummy_group = adummyvar(group);

indepvar = {'coh_prior',bsxfun(@times,coh(filt),dummy_prior(filt,:)),...
    'dot_dur_coh',bsxfun(@times,dotdur(filt),coh(filt))...
    'prior',dummy_prior(filt,:),...
    'group',dummy_group(filt,1:end-1)};
testSignificance.vars = [1,2,3,4];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

yhat = glmval(beta,x,'logit','constant','off');

% stats
fprintf(' Coh stats: %2.7f\n',LRT(1).p)
fprintf(' Prior stats: %2.7f\n',LRT(3).p)

%% Figure 2C
% Linear regression

block_prior_rel_choice = block_prior;
block_prior_rel_choice(req_choice==0) = 1-block_prior_rel_choice(req_choice==0);

% regression
filt = correct == 1;

depvar = 0.5+conf(filt)/2;
group_dummy = adummyvar(group);

indepvar = {'coh',abs(coh(filt)),...
    'prior',block_prior_rel_choice(filt),...
    'choice',choice(filt),...
    'group',group_dummy(filt,1:end)};
testSignificance.vars = [1:length(indepvar)/2];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

stats.p(idx.coh)
stats.p(idx.prior)

%% fig 4, accuracy
% logistic regression

max_t = 42;

% choice
J = trnum_eff<=max_t;
dummytr = adummyvar(trnum_eff(J));
dummygroup = adummyvar(group(J));

depvar = correct(J);
indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
    'block_prior', bsxfun(@times,block_prior_rel(J),dummytr(:,1:end)),...
    'dotdur',dotdur(J),...
    'participant',dummygroup(:,1:end-1)};

%             end
testSignificance.vars = [1,2,3,4];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

[stats.p(idx.dotdur),LRT(4).p]
LRT(3)

%% fig 4,confidence
% logistic regression

low_prc = 30;
high_conf = nan(size(conf));
for i=1:length(uni_sujs)
    K = ismember(group,uni_sujs(i));
    high_conf(K) = conf(K)>prctile(conf(K),low_prc);
end

block_prior_rel_choice = block_prior.*choice + (1-block_prior).*(1-choice);

J = trnum<=max_t;

dummytr = adummyvar(trnum(J));
dummygroup = adummyvar(group(J));

depvar = high_conf(J);

indepvar = {'coh',abs(coh(J)),'coh_group',bsxfun(@times,abs(coh(J)),dummygroup(:,1:end-1)),'trnum',dummytr,...
                'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J),...
                'coh_acc',bsxfun(@times,abs(coh(J)).*(correct(J)==0),dummygroup),...
                'dotdur',dotdur(J)};%

testSignificance.vars = [1,2,3,4,9];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

[stats.p(idx.dotdur),LRT(end).p]

%% Figure S5
% linear regression            

J = correct==1;

dummygroup = adummyvar(group);


deltab_rel_choice = d.deltab_full;
deltab_rel_choice(choice==0) = -1*deltab_rel_choice(choice==0);
belief_rel_choice = belief.*choice + (1-belief).*(1-choice);
depvar = deltab_rel_choice(J);
indepvar = {'coh',abs(coh(J)),...
                'dotdur',dotdur(J),...
                'belief',belief_rel_choice(J),...
                'participant',dummygroup(J,1:end)};
testSignificance.vars = [];
[beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar,testSignificance);

stats.p(idx.coh)
stats.p(idx.dotdur)

%%


