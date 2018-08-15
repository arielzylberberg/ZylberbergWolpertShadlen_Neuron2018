load ../data/data_CCONF.mat

fields_remove = {'choice_rel','coh_rel','block_prior_rel','trnum_back_eff','trnum'};


% for i=1:length(fields_remove)
%     dat = rmfield(dat,fields_remove{i})
% end

fields = {'choice','req_choice','coh','confidence','group','dotdur',...
    'correct','session','belief','block_prior','trnum_eff'};
renames = {'choice','req_choice','coh','confidence','group','dotdur',...
    'correct','session','belief','block_prior','trnum_eff'};

for i=1:length(dat)
    data(i).choice = dat(i).choice;
    data(i).req_choice = dat(i).req_choice;
    data(i).correct = dat(i).correct;
    data(i).coh = dat(i).coh;
    data(i).dotdur = dat(i).dotdur;
    data(i).confidence = dat(i).confidence/2 + 1/2;
    data(i).subject_number = dat(i).group;
    data(i).block_number = dat(i).session;
    data(i).belief = (dat(i).belief+1)/2;
    data(i).base_rate = dat(i).block_prior;
    data(i).trial_number = dat(i).trnum_eff;
end

save data_zylberberg_wolpert_shadlen2018 data
