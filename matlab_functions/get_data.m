function [d,m] = get_data(uni_sujs,data_path_to_file,dire_model,folder,varargin)

if nargin==3
    ff = dire_model;
elseif nargin>3
    ff = fullfile(dire_model,folder);
end

Iphi = 7;
for i=1:length(varargin)
    if isequal(varargin{i},'Iphi')
        Iphi = varargin{i+1};
    end
end

%%
nSujs = length(uni_sujs);

load(data_path_to_file);

dat = dat(uni_sujs);

for i=1:nSujs
    dat(i).belief = clip ( (dat(i).belief+1)/2 , 0 ,1);
end


if nargin>2 %asks for model data
    for i=1:nSujs
        
        m = load(fullfile(ff,['suj_',num2str(uni_sujs(i)),'.mat']));
        
        M(i).choice = m.model_allsim.choice;
        M(i).conf = m.model_allsim.conf;
        M(i).prior_belief = m.model_allsim.prior_belief;
        M(i).like = m.model_allsim.like;
        M(i).belief = m.model_allsim.belief;
        M(i).correct = bsxfun(@eq,M(i).choice,dat(i).req_choice);
        phi = m.theta(Iphi);
        
        M(i).high_conf = m.model_allsim.conf>phi;
        M(i).theta = m.theta;
        if isfield(m,'smooth') && ~isempty(m.smooth)
            f = fields(m.smooth);
            for k=1:length(f)
                M(i).smooth.(f{k}) = m.smooth.(f{k});
            end
            
        end
        %M(i).belief_scaled = ranktransfrom_one2one(M(i).belief,dat(i).belief);
        
        low_prc = m.pars.low_prc;
        dat(i).high_confidence = dat(i).confidence>prctile(dat(i).confidence,low_prc);
        
        if isfield(m,'model_fromeval')
            % M(i).like_from_conf = m.model_fromeval.right_choice_like;
            if isfield (m.model_fromeval,'right_choice_posterior_unbiased')
                M(i).like_from_conf = m.model_fromeval.right_choice_posterior_unbiased;
            end
            if isfield (m.model_fromeval,'prior_rightwards')
                M(i).prior_rightwards_fromeval = m.model_fromeval.prior_rightwards;
            end
            if isfield (m.model_fromeval,'conf_remapped')
                M(i).conf_remapped = m.model_fromeval.conf_remapped;
            end
        end
        if isfield(m,'model_fromeval') && isfield(m.model_fromeval,'K')
            nk = length(m.model_fromeval.K);
            a = bsxfun(@times,m.model_allsim.prior_belief,reshape(m.model_fromeval.K,[1,1,nk]));
            M(i).prior_rightwards_model = sum(a,3);
            M(i).K = m.model_fromeval.K;
            M(i).prior_K_ini = m.model_fromeval.prior_K_ini;
        end
        
    end
    
end

%%
clear d m

%% data

d.choice      = cat(1,dat.choice);
d.req_choice  = cat(1,dat.req_choice);
d.coh         = cat(1,dat.coh);
d.block_prior = cat(1,dat.block_prior);
d.block_prior_rel = cat(1,dat.block_prior_rel);
d.belief      = cat(1,dat.belief);
d.conf        = cat(1,dat.confidence);
if isfield(dat(1),'high_confidence')
    d.high_conf   = cat(1,dat.high_confidence);
end
d.coh_rel     = cat(1,dat.coh_rel);
d.correct     = cat(1,dat.correct);
d.trnum_back_eff = cat(1,dat.trnum_back_eff);
d.session     = dat(1).session;
for i=2:length(uni_sujs)
    d.session = [d.session; dat(i).session+nanmax(d.session)];
end
d.dotdur      = cat(1,dat.dotdur);

% comput effective trial num, ignoring nans
d.uni_session = unique(d.session);
trnum_eff = nan(size(d.session));
for i=1:length(d.uni_session)
    I = ~isnan(d.choice) & d.session==d.uni_session(i);
    trnum_eff(I) = 1:sum(I);
end
d.trnum = trnum_eff;
d.trnum_eff = trnum_eff;%duplicate, for compatibility


d.uni_block = nanunique(d.block_prior);

d.group = [];
for i=1:nSujs
    d.group = [d.group; ones(length(dat(i).confidence),1)*uni_sujs(i)];
end

d.acoh = abs(d.coh);
d.uacoh = unique(d.acoh);
d.ucoh = unique(d.coh);

d.bl = round(100*(d.block_prior.*d.req_choice+(1-d.block_prior).*(1-d.req_choice)))/100;%prior rel to req_choice

% calc delta beliefs
[d.deltab_full,d.belief_shifted] = calc_deltab(d.belief,d.session);

    
%% model
if nargin<=2
    m = nan;%return just data
else
    
    m.choice_model    = cat(1,M.choice);
    m.correct_model   = cat(1,M.correct);
    m.conf_model      =  cat(1,M.conf);
    m.conf_high_model =  cat(1,M.high_conf);
    if isfield(M(1),'like_from_conf')
        m.like_from_conf = cat(1,M.like_from_conf);% from eval
    end
    if isfield(M(1),'prior_rightwards_fromeval')
        m.prior_rightwards_fromeval = cat(1,M.prior_rightwards_fromeval);% from eval
    end
    
    if isfield(M(1),'conf_remapped')
        m.conf_remapped = cat(1,M.conf_remapped);% from eval
    end
    
    if isfield(M(1),'prior_rightwards_model')
        m.prior_rightwards_model = cat(1,M.prior_rightwards_model);
        m.K = cat(1,M.K);
        m.prior_K_ini = cat(1,M.prior_K_ini);
    end
    
    m.like_model =  cat(1,M.like);
    m.belief_model = cat(1,M.belief);
    if isfield(M(1),'right_choice_like')
        m.pright_model = cat(1,M.right_choice_like);
        pcorrect_model = pright_model;
        nn = size(pright_model,2);
        I = repmat(req_choice,1,nn)==0;
        m.pcorrect_model(I) = 1-pcorrect_model(I);
    else
        m.pright_model = m.choice_model;
        m.pcorrect_model = m.correct_model;
    end
    
    m.nsims = size(m.choice_model,2);
    
    m.smooth_flag = 0;
    if isfield(M(1),'smooth') && ~isempty(M(1).smooth)
        
        m.smooth_flag = 1;
        f = fields(M(1).smooth);
        for k=1:length(f)
            mm.(f{k}) = M(1).smooth.(f{k});
            for i=2:length(uni_sujs)
                mm.(f{k}) = cat(1,mm.(f{k}),M(i).smooth.(f{k}));
            end
        end
        mm.ucohs = M(1).smooth.ucohs;
        mm.udurs = M(1).smooth.udurs;
        
        m.ucohs = M(1).smooth.ucohs;
        m.udurs = M(1).smooth.udurs;
        
        if any(mm.ucohs==0)
            error('smooth ucohs has zeros')
        end
        
        if isfield(mm,'econf_right')
            f = {'p_right','p_right_high','p_left_high','p_right_d','p_right_high_d',...
                'p_left_high_d','econf_right','econf_left','econf_right_d','econf_left_d'};
        else
            f = {'p_right','p_right_high','p_left_high','p_right_d','p_right_high_d',...
                'p_left_high_d'};
        end
        
        for i=1:length(f)
            
            % assumes symmetric coherences
            nn = size(mm.(f{i}),2);
            i1 = 1:nn;
            i2 = (nn+1):nn;
            
            if ismatrix(mm.(f{i}))
                
                I = d.req_choice==1;
                
                mm.(f{i})(I,i1) = nan;
                I = d.req_choice==0;
                mm.(f{i})(I,i2) = nan;
                
            elseif ndims(mm.(f{i}))==3
                
                I = d.req_choice==1;
                mm.(f{i})(I,i1,:) = nan;
                I = d.req_choice==0;
                mm.(f{i})(I,i2,:) = nan;
                
            end
        end
        
        % for output
        m.smooth_model_right = mm.p_right;
        

    end
    
    
    %conf
    high_given_right = mm.p_right_high./mm.p_right;
    high_given_left = mm.p_left_high./(1-mm.p_right);
    
    [m.smooth_model_correct,...
        m.smooth_model_correct_unsigned,...
        m.smooth_model_high_given_correct,...
        m.smooth_model_high_given_error,...
        m.smooth_model_high_given_correct_unsigned,...
        m.smooth_model_high_given_error_unsigned] = ...
        sort_the_smoothies(mm.ucohs,mm.p_right,high_given_right,high_given_left,d.req_choice,d.block_prior,m.correct_model(:,1));
    
    
    m.aucohs = unique(abs(mm.ucohs));
    m.aucohs(m.aucohs<=eps) = 0; % set eps value to zero, for plotting
    
    % same, for ev
    if isfield(mm,'econf_right')
        [~,...
            ~,...
            ~,...
            ~,...
            m.smooth_model_econf_given_correct_unsigned,...
            m.smooth_model_econf_given_error_unsigned] = ...
            sort_the_smoothies(mm.ucohs,mm.p_right,mm.econf_right,mm.econf_left,d.req_choice,d.block_prior,m.correct_model(:,1));
    end
    
    
    
    if isfield(mm,'p_right_d') %multiple durations
        uc = d.uacoh;
        uc(uc==0) = 10^-9;
        uc = sort([-uc;uc]);
        if ~iseven(size(mm.p_right_high_d,2))
            error('wrong size');
        end
      
        high_given_right_d = mm.p_right_high_d./mm.p_right_d;
        high_given_left_d = mm.p_left_high_d./(1-mm.p_right_d);
        
        [m.smooth_model_correct_d,...
            m.smooth_model_correct_unsigned_d,...
            m.smooth_model_high_given_correct_d,...
            m.smooth_model_high_given_error_d,...
            m.smooth_model_high_given_correct_unsigned_d,...
            m.smooth_model_high_given_error_unsigned_d] = ...
            sort_the_smoothies( uc ,mm.p_right_d,high_given_right_d,high_given_left_d,...
            d.req_choice,d.block_prior,m.correct_model(:,1));
        
        
        % same, expected conf
        if isfield(mm,'econf_right_d')
            [~,...
                ~,...
                ~,...
                ~,...
                m.smooth_model_econf_given_correct_unsigned_d,...
                m.smooth_model_econf_given_error_unsigned_d] = ...
                sort_the_smoothies( uc ,mm.p_right_d,mm.econf_right_d,mm.econf_left_d,d.req_choice,d.block_prior,m.correct_model(:,1));
        end
        
    end
    
    %
    m.repcoh = repmat(d.coh,1,m.nsims);
    m.repdotdur = repmat(d.dotdur,1,m.nsims);
    m.repgroup = repmat(d.group,1,m.nsims);
    
    % calc delta beliefs
    m.deltab_model_full = nan(size(m.belief_model));
    m.belief_model_shifted = nan(size(m.belief_model));
    for j=1:size(m.belief_model,2)
        [m.deltab_model_full(:,j),m.belief_model_shifted(:,j)] = calc_deltab(m.belief_model(:,j),d.session);
    end
    
    
end



end

%%%%%%%%%%%%%%%%%%%%%
%%%% HELPER FN %%%%%%
%%%%%%%%%%%%%%%%%%%%%


function [smooth_model_correct,smooth_model_correct_unsigned,smooth_model_high_given_correct,smooth_model_high_given_error,...
    smooth_model_high_given_correct_unsigned,smooth_model_high_given_error_unsigned] = ...
sort_the_smoothies(ucohs,smooth_model_right,smooth_model_high_given_right,smooth_model_high_given_left,req_choice,block_prior,correct_model)

cmodel = smooth_model_right;
cmodel(:,ucohs<0,:) = 1-cmodel(:,ucohs<0,:);
if sum(ucohs==0)>0
    error('ucohs has zeros');
end

smooth_model_correct = cmodel; % out


cmodel_unsigned = cmodel;
I = req_choice==0;
cmodel_unsigned(I,:,:) = cmodel_unsigned(I,end:-1:1,:);
cmodel_unsigned(:,ucohs<0,:) = [];
smooth_model_correct_unsigned = cmodel_unsigned; %out


%conf
right_ratio = smooth_model_high_given_right;
left_ratio = smooth_model_high_given_left;

smooth_model_high_given_correct = right_ratio;
smooth_model_high_given_correct(:,ucohs<0,:) = left_ratio(:,ucohs<0,:);

I = correct_model == 0;
smooth_model_high_given_correct(I,:,:) = nan; % out

smooth_model_high_given_error = right_ratio;
smooth_model_high_given_error(:,ucohs>0,:) = left_ratio(:,ucohs>0,:);

I = correct_model == 1;
smooth_model_high_given_error(I,:,:) = nan; %out

aux = smooth_model_high_given_correct;
aux(req_choice==0,:,:) = aux(req_choice==0,end:-1:1,:);
smooth_model_high_given_correct_unsigned = aux(:,ucohs>=0,:);

aux = smooth_model_high_given_error;
aux(req_choice==0,:,:) = aux(req_choice==0,end:-1:1,:);
smooth_model_high_given_error_unsigned = aux(:,ucohs>=0,:);


end