classdef plot_functions < handle
    properties
        
        d
        m
        
        uni_sujs
        nsujs
        
        G
        titles
        colores_us
        colores_uscoh
        
        figures_handle
        p_handle
        
        figpos % figure positions
        
        
        
    end
    
    methods
        
        
        function obj = plot_functions(d,m)
            
            obj.d = d;
            obj.m = m;
            
            uni_sujs = nanunique(d.group);
            obj.G{1} = uni_sujs;
            obj.titles{1} = 'all_sujs';
            for i=1:length(uni_sujs)
                obj.G{i+1} = uni_sujs(i);
                obj.titles{i+1} = ['S',num2str(uni_sujs(i))];
            end
            
            aux_colors = cbrewer('qual','Paired',length(d.uni_block));
            obj.colores_us = aux_colors([1,3,5,6,4,2],:);
            
            obj.colores_uscoh = rainbow_colors(6,'colorType',6);
            
            obj.uni_sujs = uni_sujs;
            obj.nsujs = length(uni_sujs);
            
            obj.d.block_bias = max(obj.d.block_prior,1-obj.d.block_prior);
            color_block_bias = movshon_colors(length(unique(obj.d.block_bias)));
            obj.m.color_block_bias = color_block_bias;
            
            uni_prior_rel_choice = unique(obj.d.bl);
            obj.d.uni_prior_rel_choice = uni_prior_rel_choice;
            
            
            obj.m.coloresm = movshon_colors(length(uni_prior_rel_choice));
            
            
            obj.figpos.onen = [119   351  1715*(length(uni_sujs)+1)/6   258];
            obj.figpos.twon = [119   351  1715*(length(uni_sujs)+1)/6   500];
            
            obj.figpos.onen_u = [1   578  1440*(length(uni_sujs)+1)/6   184];
            obj.figpos.twon_u = [1   578  1440*(length(uni_sujs)+1)/6   368];
            
            % remap
            obj.d.unscaled_belief = d.belief;
            obj.d.unscaled_delta_belief = d.deltab_full;
            
        end
        
        
        function p = likelihood_distributions_per_coh(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            nuacoh = length(uacoh);
            uni_group = unique(group);
            
            sujs_to_include = uni_group(1);
            nn = length(sujs_to_include);
            
            p = publish_plot(nuacoh,nn);
            set(gcf,'position',[570   175   430   623]);
            colores = rainbow_colors(nuacoh);
            
            like_model_R = like_model;
            like_model_R(choice_model==0) = 1-like_model_R(choice_model==0);
            for i=1:nuacoh
                for j=1:nn
                    I = acoh == uacoh(i) & group==sujs_to_include(j);
                    p.next();
                    edges = linspace(0,1,100);
                    h = histc(like_model_R(I,1),edges);
                    hb = bar(edges,h,'histc');
                    set(hb,'facecolor',colores(i,:));
                    ht(i,j) = title(['coh=',num2str(uacoh(i)*100),'%']);
                end
            end
            p.format('FontSize',14);
            xlabel('Counterfactual confidence (p_u(R))')
            set(p.h_ax,'xlim',[0,1],'tickdir','out','ycolor','none')
            p.unlabel_center_plots();
            set(ht,'FontWeight','normal');
            
        end
        
        function p = dbelief_vs_like_from_model(obj,varargin)
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = length(obj.uni_sujs);
            
            belief_range = [0,1];
            for i=1:length(varargin)
                if isequal(varargin{i},'belief_range')
                    belief_range = varargin{i+1};
                end
            end
            
            p = publish_plot(2,nSujs+1);
            
            moving_av_flag = 1;
            if moving_av_flag
                
                ntr_av = 100;
                
                % data
                for ig=1:(nSujs+1)
                    I = ismember(obj.d.group,obj.G{ig}) & in_range(belief,belief_range);
                    [X,Y] = dotsanalysis.choice_against_duration_moving(deltab_full,choice,like_from_conf,ntr_av, 0 ,I);
                    
                    p.next();
                    plot(X{1},Y{1},'linewidth',1);
                    hold all
                    plot(X{2},Y{2},'linewidth',1);
                    title(obj.titles{ig},'interpreter','none')
                    if ig==1
                        ylabel('\Delta belief,data')
                    end
                end
                
                % model
                add_noise_to_like_flag = 0; % adds gaussian noise and clips
                if add_noise_to_like_flag
                    like_model = clip(like_model + randn(size(like_model))*0.1,0,1);
                end
                m_like_right = like_model;
                m_like_right(choice_model==0) = 1-m_like_right(choice_model==0);
                
                % recalc like - doesnt work, need prior_K for the model
                %                 [m_like_right,~] = ...
                %                     calc_like_from_conf(prior_K, K, choice_model(:,1), conf_model(:,1));
                
                for ig=1:(nSujs+1)
                    I = ismember(obj.d.group,obj.G{ig}) & in_range(belief_model(:,1),belief_range);
                    p.next();
                    [X,Y] = dotsanalysis.choice_against_duration_moving(deltab_model_full(:,1),...
                        choice_model(:,1),m_like_right,ntr_av, 0 ,I);
                    
                    plot(X{1},Y{1},'linewidth',1);
                    hold all
                    plot(X{2},Y{2},'linewidth',1);
                    
                    title(obj.titles{ig},'interpreter','none')
                    if ig==1
                        ylabel('\Delta belief,model')
                    end
                    xlabel('counterfactual conf')
                end
                
            else
                
                % data
                [~,idx,~,idx_val] = index_prctile(like_from_conf+randn(size(choice))*0.0000001,[0:5:100]);
                for ig=1:(nSujs+1)
                    I = ismember(obj.d.group,obj.G{ig}) & in_range(belief,belief_range);
                    p.next();
                    for iChoice=0:1
                        curva_media(deltab_full,idx_val,choice==iChoice & I,2);
                        hold all
                    end
                    title(obj.titles{ig},'interpreter','none')
                    if ig==1
                        ylabel('\Delta belief,data')
                    end
                end
                
                % model
                add_noise_to_like_flag = 0; % adds gaussian noise and clips
                if add_noise_to_like_flag
                    like_model = clip(like_model + randn(size(like_model))*0.1,0,1);
                end
                m_like_right = like_model;
                m_like_right(choice_model==0) = 1-m_like_right(choice_model==0);
                [~,idx,~,idx_val] = index_prctile(m_like_right(:,1)+randn(size(choice))*0.0000001,[0:5:100]);
                for ig=1:(nSujs+1)
                    I = ismember(obj.d.group,obj.G{ig}) & in_range(belief_model(:,1),belief_range);
                    p.next();
                    for iChoice=0:1
                        curva_media(deltab_model_full(:,1),idx_val,choice_model(:,1)==iChoice & I,2);
                        hold all
                    end
                    title(obj.titles{ig},'interpreter','none')
                    if ig==1
                        ylabel('\Delta belief,model')
                    end
                    xlabel('counterfactual conf')
                end
            end
            
            p.unlabel_center_plots();
            set(p.h_ax,'xlim',[0,1]);
            same_ylim(p.h_ax);
            for i=1:length(p.h_ax)
                p.current_ax(i);
                hold all
                plot(xlim,[0,0],'color',0.7*[1,1,1],'linestyle',':')
                plot([0.5,0.5],ylim,'color',0.7*[1,1,1],'linestyle',':')
            end
            
            %set(p.h_ax,'xlim',[0,1])
            p.format();
            set(gcf,'Position',obj.figpos.twon)
            %nontouching_spines(p.h_ax);
            
            % copy the all suj to new fig
            %             pp = publish_plot(1,2);
            %             set(gcf,'Position',[342  350  582  289]);
            %             pp.copy_from_ax(p.h_ax(1),1);
            %             pp.copy_from_ax(p.h_ax(5),2);
            %             pp.save('ancho',12);
            
            
        end
        
        function p = regression_choice_conf(obj,nsimu)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            if nargin==1 || isempty(nsimu)
                nsimu = nsims;
            end
            
            p = publish_plot(2,1);
            p.next();
            
            max_t = 42;
            I = ismember(group,obj.G{ 1 });
            J = I & trnum_eff<=max_t;
            
            dummytr = adummyvar(trnum_eff(J));
            dummygroup = adummyvar(group(J));
            
            depvar = correct(J);
            indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
                'block_prior', bsxfun(@times,block_prior_rel(J),dummytr(:,1:end)),...
                'dotdur',dotdur(J),... %added
                'participant',dummygroup(:,1:end-1)};
            
            
            
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
            ww_data = beta(idx.block_prior);
            se_data = stats.se(idx.block_prior);
            
            ww_model = nan(max_t,nsims);
            se_model = nan(max_t,nsims);
            nreg = nsims;
            for k=1:nreg % different simulations
                
                depvar = correct_model(J,k);
                
                [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
                
                ww_model(:,k) = beta(idx.block_prior);
                se_model(:,k) = stats.se(idx.block_prior);
            end
            
            calc_errors_from_sims = 1;
            if calc_errors_from_sims
                K = abs(zscore(ww_model(:)))>20;%some didnt converge
                ww_model(K) = nan;
                e_model = nanstd(ww_model,[],2);
            else
                e_model = nanmean(se_model,2) * 1.96;
            end
            
            a = [nanmean(ww_model,2)-e_model, nanmean(ww_model,2)+e_model];
            xx = [1:max_t,max_t:-1:1];
            a(:,2) = a(end:-1:1,2);
            patch(xx,a(:),0.8*[1,1,1],'edgecolor','none');
            hold all
            %terrorbar(1:length(ww_data),ww_data,se_data,'LineStyle','none','color','k','marker','o','markersize',8,'markerfacecolor',0.7*[1,1,1])
            terrorbar(1:max_t,ww_data,se_data,'LineStyle','none','color','k','marker','none')
            plot(1:max_t,ww_data,'LineStyle','none','color','k','marker','o','markersize',6,'markerfacecolor',0.7*[1,1,1])
            
            %             terrorbar(1:length(ww_data),ww_data,se_data,'LineStyle','-','color','k','marker','o','markersize',4,'markerfacecolor','k')
            %             hold all
            %             plot(1:size(ww_model,1),nanmean(ww_model,2),'LineStyle','--','color','k','marker','o','markersize',4,'markerfacecolor','w')
            plot([1,max_t],[0,0],'color',0*[1,1,1])
            
            %             p.format();
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_betas_bias_influence'));
            
            %             xlabel('Trial number')
            ylabel('\beta, choice')
            %             obj.figures_handle.regression_choice = p.h_fig;
            %pchoice = p;
            
            
            
            %end
            
            
            %function regression_conf(obj)
            %             struct2vars(obj.d);
            %             struct2vars(obj.m);
            
            block_prior_rel_choice = block_prior.*choice+(1-block_prior).*(1-choice);
            block_prior_rel_choice_model = bsxfun(@times,block_prior,choice_model)+...
                bsxfun(@times,1-block_prior,1-choice_model);
            
            %p = publish_plot(1,1);
            
            I = ismember(group,obj.G{ 1 });
            
            max_t = 42;
            
            
            p.next();
            
            
            %depvar = conf(J)/2 + 0.5;
            
            %             indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
            %                 'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
            %                 'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J)};%,...
            
            %             indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
            %                 'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
            %                 'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J),...
            %                 'dotdur',dotdur(J)};
            
            
            conf_regre = 1;% initial neuron submission
            
            %             if exist('conf_remapped','var')
            %                 % use analog conf
            %                 conf_regre = 5;
            %             end
            
            switch conf_regre
                case 1
                    J = I & trnum<=max_t;
                    depvar = high_conf(J);
                    dummytr = adummyvar(trnum(J));
                    dummygroup = adummyvar(group(J));
                    indepvar = {'coh',abs(coh(J)),'coh_group',bsxfun(@times,abs(coh(J)),dummygroup(:,1:end-1)),'trnum',dummytr,...
                        'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                        'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J),...
                        'coh_x_acc',bsxfun(@times,abs(coh(J)).*(correct(J)==1),dummygroup),...
                        'dotdur',dotdur(J)};%
                case 2
                    J = I & trnum<=max_t;
                    depvar = high_conf(J);
                    dummytr = adummyvar(trnum(J));
                    dummygroup = adummyvar(group(J));
                    indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
                        'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                        'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J)};%,...
                case 3
                    J = I & trnum<=max_t;
                    depvar = high_conf(J);
                    dummytr = adummyvar(trnum(J));
                    dummygroup = adummyvar(group(J));
                    indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
                        'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                        'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J),'dotdur',dotdur(J)};%,...
                    
                case 4
                    J = I & trnum<=max_t & correct==1;
                    depvar = high_conf(J);
                    dummytr = adummyvar(trnum(J));
                    dummygroup = adummyvar(group(J));
                    indepvar = {'coh',bsxfun(@times,abs(coh(J)),dummygroup),'trnum',dummytr,...
                        'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                        'participant',dummygroup(:,1:end-1),'choice',choice(J),'dotdur',dotdur(J)};%,...
                    
                    
                case 5
                    J = I & trnum<=max_t;
                    depvar = conf_remapped(J);
                    dummytr = adummyvar(trnum(J));
                    dummygroup = adummyvar(group(J));
                    indepvar = {'coh',abs(coh(J)),'coh_group',bsxfun(@times,abs(coh(J)),dummygroup(:,1:end-1)),'trnum',dummytr,...
                        'block_prior', bsxfun(@times,block_prior_rel_choice(J),dummytr),...
                        'participant',dummygroup(:,1:end-1),'choice',choice(J),'correct',correct(J),...
                        'coh_x_acc',bsxfun(@times,abs(coh(J)).*(correct(J)==1),dummygroup),...
                        'dotdur',dotdur(J)};%
                    
            end
            
            
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
            ww_data = beta(idx.block_prior);
            se_data = stats.se(idx.block_prior);
            
            % model
            
            ww_model = nan(max_t,nsims);
            se_model = nan(max_t,nsims);
            
            for k=1:nsimu % different simulations
                
                
                
                %depvar = conf_model(J,k);
                
                
                %before
                
                
                
                J = I & trnum<=max_t;
                depvar = conf_high_model(J,k);
                dummytr = adummyvar(trnum(J));
                dummygroup = adummyvar(group(J));
                indepvar = {'coh',abs(coh(J)),'coh_group',bsxfun(@times,abs(coh(J)),dummygroup(:,1:end-1)),'trnum',dummytr,...
                    'block_prior', bsxfun(@times,block_prior_rel_choice_model(J,k),dummytr),...
                    'participant',dummygroup(:,1:end-1),'choice',choice_model(J,k),'correct',correct_model(J,k),...
                    'coh_x_acc',bsxfun(@times,abs(coh(J)).*(correct_model(J,k)==1),dummygroup),...
                    'dotdur',dotdur(J)};%
                
                
                [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
                ww_model(:,k) = beta(idx.block_prior);
                se_model(:,k) = stats.se(idx.block_prior);
            end
            
            if calc_errors_from_sims
                e_model = nanstd(ww_model,[],2);
            else
                e_model = nanmean(se_model,2) * 1.96;
            end
            
            a = [nanmean(ww_model,2)-e_model, nanmean(ww_model,2)+e_model];
            xx = [1:max_t,max_t:-1:1];
            a(:,2)=a(end:-1:1,2);
            patch(xx,a(:),0.8*[1,1,1],'edgecolor','none');
            hold all
            %plot(1:max_t,nanmean(ww_model,2),'r')
            terrorbar(1:max_t,ww_data,se_data,'LineStyle','none','color','k','marker','none')
            plot(1:max_t,ww_data,'LineStyle','none','color','k','marker','o','markersize',6,'markerfacecolor',0.7*[1,1,1])
            %plot(1:size(ww_model,1),nanmean(ww_model,2),'LineStyle','--','color','k','marker','o','markersize',4,'markerfacecolor','w')
            hold all
            plot([1,max_t],[0,0],'color',0*[1,1,1])
            xlabel('Trial number')
            ylabel('\beta, confidence')
            p.format();
            
            set(p.h_ax,'xlim',[0,43]);
            set(p.h_ax(1),'xticklabel','');
            p.displace_ax(2,0.1,2);
            set(gcf,'Position',[507  145  493  653])
            
            
            set(p.h_ax(1),'ylim',[-2,10],'xcolor','none');
            set(p.h_ax(2),'ylim',[-1,5]);
            
            % nontouching_spines(p.h_ax);
            
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_betas_bias_influence_conf'));
            obj.figures_handle.regression_conf = p.h_fig;
            %pconf = p;
            
            
        end
        
        
        function p = combined_fig4_opt2(obj)
            
            use_analog_confidence = 0;
            
            p = publish_plot(2,3);
            set(gcf,'Position',[187  264  854  423])
            p.resize_horizontal(1:3,[1.2,1,1.5])
            p.resize_horizontal(4:6,[1.2,1,1.5])
            
            p.displace_ax([4:6],.08,2);
            p.displace_ax([1,4],.03,1);
            
            p.displace_ax([3,6],.05,1);
            
            paux = obj.acc_unsigned();
            p.copy_from_ax(paux.h_ax(1),1,1);
            
            if use_analog_confidence
                paux = obj.conf_continuous_unsigned_per_suj();
                p.copy_from_ax(paux.h_ax(1),4,1);
            else
                paux = obj.conf_unsigned_per_suj();
                p.copy_from_ax(paux.h_ax(1),4,1);
            end
            
            paux = obj.accuracy_vs_duration();
            p.copy_from_ax(paux.h_ax(1),2,1);
            
            if use_analog_confidence
                paux = obj.confidence_analog_correct_vs_duration_split_by_coh();
                p.copy_from_ax(paux.h_ax(1),5,1);
            else
                paux = obj.confidence_correct_vs_duration_split_by_coh();
                p.copy_from_ax(paux.h_ax(1),5,1);
            end
            
            BreakXAxis(p.h_ax([1,4]),0.18,'break_scaling',4);
            
            paux = obj.regression_choice_conf();
            p.copy_from_ax(paux.h_ax(1),3, 0);
            p.copy_from_ax(paux.h_ax(2),6, 1);
            
            set(p.h_ax(1:3),'xticklabel','');
            for i=1:3
                set(get(p.h_ax(i),'xlabel'),'String','');
            end
            set(p.h_ax,'XMinorTick','off','YMinorTick','off');
            set(p.h_ax(1:2),'ylim',[.4,1],'ytick',0.4:0.2:1);
            if use_analog_confidence
                set(p.h_ax(4:5),'ylim',[.6,1],'ytick',0.6:0.2:1);
            else
                set(p.h_ax(4:5),'ylim',[.2,1],'ytick',0.2:0.4:1);
            end
            set(p.h_ax([2,5]),'ycolor',0.7*[1,1,1])
            set(get(p.h_ax(2),'ylabel'),'String','');
            set(get(p.h_ax(5),'ylabel'),'String','');
            set(p.h_ax([2,5]),'yticklabel','');
            
            p.current_ax(1);
            ylabel('Proportion correct');
            
            p.current_ax(5);
            xlabel('Motion duration (s)')
            
            %p.current_ax(2);
            %ylabel('Proportion correct');
            
            p.current_ax(4);
            xlabel('Motion strength (%coh)')
            
            if use_analog_confidence
                ylabel('Confidence');
            else
                ylabel({'Proportion','high confidence'});
            end
            
            
            %p.current_ax(5);
            %ylabel('P high confidence');
            
            
            p.current_ax(3);
            ylabel({'Base-rate leverage','on accuracy'})
            %title('Accuracy')
            
            p.current_ax(6);
            ylabel({'Base-rate leverage','on confidence'})
            %title('Confidence')
            
            
            
            h = findobj(p.h_fig,'Type','line');
            set(h,'LineWidth',.7);
            
            drawnow
            
        end
        
        
        
        
        function choice_vs_belief_equispaced(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            %perf_include = [1];
            perf_include = [0,1];
            %perf_include = [0];
            
            deltab = choice;
            deltab_model = choice_model;
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            
            p = publish_plot(2,1);
            set(gcf,'Position',[928   109   312   657])
            p.shrink([1,2],.8,.8)
            %p.displace_ax(2,0.08,2);
            
            
            nbins = 10;
            edges = linspace(0,1,nbins);
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            bins_data = discretize(belief_shifted,edges);
            tt = edges(1:end-1)+diff(edges)/2;
            colores = rainbow_colors(length(ucoh));
            colores(ucoh==0,:) = [0,0,0];
            for i=1:length(ucoh)
                
                I = coh==ucoh(i) & ~isnan(deltab) & S;
                if sum(I)>0
                    [~,xx,ss] = curva_media(deltab,bins_data,I,0);
                    hold all
                    hp(i) = terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle','-','linewidth',lw);
                end
                
            end
            
            ylabel('P rightward')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(2);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            bins_model = discretize(belief_model_shifted(:),edges);
            
                for i=1:length(ucoh)
                    I = repcoh(:)==ucoh(i) & ~isnan(deltab_model(:)) & SS;
                    %[tt,xx,ss] = dotsanalysis.choice_against_duration(deltab_model(:),ones(nn,1),belief_model_shifted(:),nbins,0,I);
                    [~,xx,ss] = curva_media(deltab_model(:),bins_model,I,0);
                    hold all
                    terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle','-','linewidth',lw);
                end
            
            
            xlabel('Current belief')
            ylabel('P rightward')
            
            
            %npoints = npoints_global*3;
            
           
            %ht(2) = title('Model');
            
            for i=1:2
                p.current_ax(i);
                hold all
                plot([0,1],[0.5,0.5],'k--')
            end
            
            set(p.h_ax,'ylim',[0,1],'ytick',[0:0.25:1],'xlim',[0,1])
            
            p.current_ax(1);
            
            
            set(p.h_ax,'tickdir','out','color','none','xtick',[0:0.25:1]);
            %set(ht,'fontweight','normal')
            
            p.unlabel_center_plots();
            
            p.format('FontSize',16);
            
            hl = legend_n(ucoh*100);
            set(hl,'Location','NorthEastOutside','FontSize',10);
            hll = get(hl,'title');
            set(hll,'string','coh (%)');
            
            drawnow
            p.displace_ax(2,0.08,2);
            
            p.current_ax(1);
            title('Data')
            p.current_ax(2);
            title('Bayesian model')
            
        end
        
        function dbelief_vs_belief_equispaced(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            perf_include = [1];
            % perf_include = [0,1];
            %perf_include = [0];
            
            deltab = deltab_full;
            deltab_model = deltab_model_full;
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            
            p = publish_plot(2,2);
            set(gcf,'Position',[1036   845   524   493])
            p.displace_ax([3,4],0.08,2);
            
            
            nbins = 10;
            edges = linspace(0,1,nbins);
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            bins_data = discretize(belief_shifted,edges);
            tt = edges(1:end-1)+diff(edges)/2;
            for i=1:length(uacoh)
                for k=1:2
                    I = acoh==uacoh(i) & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [~,xx,ss] = curva_media(deltab,bins_data,I,0);
                        hold all
                        hp(i) = terrorbar(tt,xx,ss,'color',obj.colores_uscoh(i,:),'LineStyle',lsty{k},'linewidth',lw);
                    end
                end
            end
            
            ylabel('Change in belief')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(3);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            bins_model = discretize(belief_model_shifted(:),edges);
            for k=1:2
                for i=1:length(uacoh)
                    I = abs(repcoh(:))==uacoh(i) & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    %[tt,xx,ss] = dotsanalysis.choice_against_duration(deltab_model(:),ones(nn,1),belief_model_shifted(:),nbins,0,I);
                    [~,xx,ss] = curva_media(deltab_model(:),bins_model,I,0);
                    hold all
                    terrorbar(tt,xx,ss,'color',obj.colores_uscoh(i,:),'LineStyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief')
            
            
            %npoints = npoints_global*3;
            
            p.current_ax(2);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            
            [uni_dur,Idur] = index_prctile(dotdur,[0,50,100]);
            % colores = [0,0,1;1,0,0];
            colores = [0.1529,0.6667,0.8824;1,0,0];
            for i=1:length(uni_dur)-1
                for k=1:2
                    I = Idur==i & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [~,xx,ss] = curva_media(deltab,bins_data,I,0);
                        hold all
                        hp(i) = terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                        
                    end
                end
            end
            
            ylabel('Change in belief; Data')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(4);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            Idurrep = repmat(Idur,nsims,1);
            for k=1:2
                for i=1:length(uni_dur)-1
                    I = Idurrep==i & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    %[tt,xx,ss] = dotsanalysis.choice_against_duration(deltab_model(:),ones(nn,1),belief_model_shifted(:),nbins,0,I);
                    hold all
                    %terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                    [~,xx,ss] = curva_media(deltab_model(:),bins_model,I,0);
                    hp(i) = terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief; Model')
            %ht(2) = title('Model');
            
            for i=1:4
                p.current_ax(i);
                hold all
                plot([0,1],[0,0],'k--')
            end
            
            set(p.h_ax,'ylim',[-0.15,0.15],'ytick',[-0.1,0,0.1],'xlim',[0,1])
            
            p.current_ax(1);
            
            
            set(p.h_ax,'tickdir','out','color','none','xtick',[0:0.25:1]);
            %set(ht,'fontweight','normal')
            
            p.unlabel_center_plots();
            
            p.format('FontSize',16);
            
            drawnow
            
        end
        
        
        function dbelief_vs_belief_not_moving_av(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            perf_include = [1];
            % perf_include = [0,1];
            %perf_include = [0];
            
            deltab = deltab_full;
            deltab_model = deltab_model_full;
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            
            p = publish_plot(2,2);
            set(gcf,'Position',[1036   845   524   493])
            p.displace_ax([3,4],0.08,2);
            
%             if perf_include==0
%                 npoints_global = 20;
%             else
%                 npoints_global = 200;
%             end
%             npoints = npoints_global;
            
            nbins = 20;
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            noise_for_prc = 10^-10*rand(size(belief_shifted)); % hack to compute percentiles
            for i=1:length(uacoh)
                for k=1:2
                    I = acoh==uacoh(i) & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [tt,xx,ss] = dotsanalysis.choice_against_duration(deltab,ones(size(coh)),belief_shifted + ...
                            noise_for_prc,nbins,0,I);
                        hold all
                        hp(i) = terrorbar(tt,xx,ss,'color',obj.colores_uscoh(i,:),'LineStyle',lsty{k},'linewidth',lw);
                        
                    end
                end
            end
            
            ylabel('Change in belief')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(3);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            for k=1:2
                for i=1:length(uacoh)
                    I = abs(repcoh(:))==uacoh(i) & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    [tt,xx,ss] = dotsanalysis.choice_against_duration(deltab_model(:),ones(nn,1),belief_model_shifted(:),nbins,0,I);
                    hold all
                    terrorbar(tt,xx,ss,'color',obj.colores_uscoh(i,:),'LineStyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief')
            
            
            %npoints = npoints_global*3;
            
            p.current_ax(2);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            
            [uni_dur,Idur] = index_prctile(dotdur,[0,50,100]);
            % colores = [0,0,1;1,0,0];
            colores = [0.1529,0.6667,0.8824;1,0,0];
            noise_for_prc = 10^-10*rand(size(belief_shifted)); % hack to compute percentiles
            for i=1:length(uni_dur)-1
                for k=1:2
                    I = Idur==i & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [tt,xx,ss] = dotsanalysis.choice_against_duration(deltab,ones(size(coh)),belief_shifted + ...
                            noise_for_prc,nbins,0,I);
                        hold all
                        terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                        
                    end
                end
            end
            
            ylabel('Change in belief; Data')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(4);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            Idurrep = repmat(Idur,nsims,1);
            for k=1:2
                for i=1:length(uni_dur)-1
                    I = Idurrep==i & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    [tt,xx,ss] = dotsanalysis.choice_against_duration(deltab_model(:),ones(nn,1),belief_model_shifted(:),nbins,0,I);
                    hold all
                    terrorbar(tt,xx,ss,'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief; Model')
            %ht(2) = title('Model');
            
            for i=1:4
                p.current_ax(i);
                hold all
                plot([0,1],[0,0],'k--')
            end
            
            set(p.h_ax,'ylim',[-0.15,0.15],'ytick',[-0.1,0,0.1],'xlim',[0,1])
            
            p.current_ax(1);
            
            
            set(p.h_ax,'tickdir','out','color','none','xtick',[0:0.25:1]);
            %set(ht,'fontweight','normal')
            
            p.unlabel_center_plots();
            
            p.format('FontSize',16);
            
            drawnow
            
        end
        
        
        %keep
        function dbelief_vs_belief(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            perf_include = [1];
            % perf_include = [0,1];
            %perf_include = [0];
            
            deltab = deltab_full;
            deltab_model = deltab_model_full;
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            
            p = publish_plot(2,2);
            set(gcf,'Position',[1036   845   524   493])
            p.displace_ax([3,4],0.08,2);
            
            if perf_include==0
                npoints_global = 20;
            else
                npoints_global = 200;
            end
            
            npoints = npoints_global;
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            for i=1:length(uacoh)
                for k=1:2
                    I = acoh==uacoh(i) & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),belief_shifted,npoints,0,I);
                        hold all
                        if ~isempty(tt{1})
                            hp(i) = plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:),'LineStyle',lsty{k},'linewidth',lw);
                        end
                    end
                end
            end
            
            ylabel('Change in belief')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(3);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            for k=1:2
                for i=1:length(uacoh)
                    I = abs(repcoh(:))==uacoh(i) & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),belief_model_shifted(:),npoints * nsims,0,I);
                    hold all
                    plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:),'linestyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief')
            
            
            npoints = npoints_global*3;
            
            p.current_ax(2);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            lsty = {'--','-'};
            
            [uni_dur,Idur] = index_prctile(dotdur,[0,50,100]);
            % colores = [0,0,1;1,0,0];
            colores = [0.1529,0.6667,0.8824;1,0,0];
            for i=1:length(uni_dur)-1
                for k=1:2
                    I = Idur==i & ~isnan(deltab) & S & req_choice==uni_req(k);
                    if sum(I)>0
                        [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),belief_shifted,npoints,0,I);
                        hold all
                        if ~isempty(tt{1})
                            hp(i) = plot(tt{1},xx{1},'color',colores(i,:),'LineStyle',lsty{k},'linewidth',lw);
                        end
                    end
                end
            end
            
            ylabel('Change in belief; Data')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(4);
            
            SS = repmat(S,nsims,1);
            rep_req_choice = repmat(req_choice,nsims,1);
            nn = length(repcoh(:));
            Idurrep = repmat(Idur,nsims,1);
            for k=1:2
                for i=1:length(uni_dur)-1
                    I = Idurrep==i & ~isnan(deltab_model(:)) & SS & rep_req_choice == uni_req(k);
                    [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),belief_model_shifted(:),npoints * nsims,0,I);
                    hold all
                    plot(tt{1},xx{1},'color',colores(i,:),'linestyle',lsty{k},'linewidth',lw);
                end
            end
            
            xlabel('Current belief')
            ylabel('Change in belief; Model')
            %ht(2) = title('Model');
            
            for i=1:4
                p.current_ax(i);
                hold all
                plot([0,1],[0,0],'k--')
            end
            
            set(p.h_ax,'ylim',[-0.15,0.15],'ytick',[-0.1,0,0.1])
            
            p.current_ax(1);
            
            
            set(p.h_ax,'tickdir','out','color','none','xtick',[0:0.25:1]);
            %set(ht,'fontweight','normal')
            
            p.unlabel_center_plots();
            
            p.format('FontSize',16);
            
            drawnow
            
        end
        
        
        % for R2, Neuron
        function choice_vs_belief(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            % perf_include = [1];
            perf_include = [0,1];
            %perf_include = [0];
            
            deltab = choice; % notation hack
            deltab_model = choice_model; %notation hack
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            
            p = publish_plot(2,2);
            set(gcf,'Position',[1036   845   524   493])
            p.displace_ax([3,4],0.08,2);
            
            if perf_include==0
                npoints_global = 20;
            else
                npoints_global = 200;
            end
            
            npoints = npoints_global;
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            for i=1:length(uacoh)
                I = acoh==uacoh(i) & ~isnan(deltab) & S;
                if sum(I)>0
                    [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),belief_shifted,npoints,0,I);
                    hold all
                    if ~isempty(tt{1})
                        hp(i) = plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:),'LineStyle','-','linewidth',lw);
                    end
                end
            end
            
            ylabel('p rightward')
            %ht(1) = title('Data');
            
            
            
            % model
            p.current_ax(3);
            
            SS = repmat(S,nsims,1);
            nn = length(repcoh(:));
                for i=1:length(uacoh)
                    I = abs(repcoh(:))==uacoh(i) & ~isnan(deltab_model(:)) & SS;
                    [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),belief_model_shifted(:),npoints * nsims,0,I);
                    hold all
                    plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:),'linestyle','-','linewidth',lw);
                end
            
            
            xlabel('Current belief')
            ylabel('p rightward')
            
            
            npoints = npoints_global*3;
            
            p.current_ax(2);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            uni_req = [0,1];
            %lsty = {'--','-'};
            
            [uni_dur,Idur] = index_prctile(dotdur,[0,50,100]);
            % colores = [0,0,1;1,0,0];
            colores = [0.1529,0.6667,0.8824;1,0,0];
            for i=1:length(uni_dur)-1
                
                    I = Idur==i & ~isnan(deltab) & S;
                    if sum(I)>0
                        [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),belief_shifted,npoints,0,I);
                        hold all
                        if ~isempty(tt{1})
                            hp(i) = plot(tt{1},xx{1},'color',colores(i,:),'LineStyle','-','linewidth',lw);
                        end
                    end
                
            end
            
            ylabel('p rightward; Data')
            %ht(1) = title('Data');
            
            % model
            p.current_ax(4);
            
            SS = repmat(S,nsims,1);
            nn = length(repcoh(:));
            Idurrep = repmat(Idur,nsims,1);
            
            for i=1:length(uni_dur)-1
                I = Idurrep==i & ~isnan(deltab_model(:)) & SS;
                [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),belief_model_shifted(:),npoints * nsims,0,I);
                hold all
                plot(tt{1},xx{1},'color',colores(i,:),'linestyle','-','linewidth',lw);
            end
            
            xlabel('Current belief')
            ylabel('p rightward; Model')
            %ht(2) = title('Model');
            
            for i=1:4
                p.current_ax(i);
                hold all
                plot([0,1],[0.5,0.5],'k--')
            end
            
            set(p.h_ax,'ylim',[0,1],'ytick',[0,1])
            
            p.current_ax(1);
            
            
            set(p.h_ax,'tickdir','out','color','none','xtick',[0:0.25:1]);
            %set(ht,'fontweight','normal')
            
            p.unlabel_center_plots();
            
            p.format('FontSize',16);
            
            drawnow
            
        end
        
        %keep
        function p = dbelief_vs_belief_unsplit_by_choice_split_prev_coh(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            %perf_include = [1];
            perf_include = [0,1];
            %perf_include = [0];
            
            deltab = deltab_full;
            deltab_model = deltab_model_full;
            
            % filter
            deltab(~ismember(correct,perf_include)) = nan;
            deltab_model(~ismember(correct_model,perf_include)) = nan; % only correct trials
            
            p = publish_plot(1,2);
            set(gcf,'Position',[254  365  707  300])
            
            
            npoints_global = 20;
            %                 if perf_include==0
            %                     npoints_global = 20;
            %                 else
            %                     npoints_global = 200;
            %
            %                     %npoints_global = 60;
            %                 end
            
            npoints = npoints_global;
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            lw = 1;
            S = ismember(group,obj.G{ig});
            
            
            coh_prev = nan(size(coh));
            coh_prev(2:end) = coh(1:end-1);
            coh_prev(trnum_eff==1) = nan;
            
            lsty = '-';
            % colores = [0,0,1;0,0,0;1,0,0];
            colores = [0,0,1;1,0,0];
            [~,~,~,v] = index_prctile(belief_shifted + randn(size(belief_shifted))*0.00001,[0:20:100]);
            %cats = [coh_prev<0,coh_prev==0,coh_prev>0];
            cats = [coh_prev<0,coh_prev>0];
            %cats = [coh_prev==-0.256,coh_prev==0,coh_prev==0.256];
            for i=1:size(cats,2)
                %for k=1:2
                I = coh==0 & cats(:,i)==1 & ~isnan(deltab) & S;
                % I = coh==0 & cats(:,i)==1 & ~isnan(deltab) & S & trnum_eff==2;
                %I = acoh==uacoh(i) & ~isnan(deltab) & S & req_choice==uni_req(k);
                if sum(I)>0
                    [tt,xx,ss] = curva_media(deltab,v,I,0);
                    terrorbar(tt,xx,ss,'marker','.','linestyle','-','color',colores(i,:),'linewidth',lw);
                    hold all
                    %                             [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),belief_shifted,npoints,0,I);
                    %                             hold all
                    %                             if ~isempty(tt{1})
                    %                                 hp(i) = plot(tt{1},xx{1},'color',colores(i,:),'LineStyle',lsty,'linewidth',lw);
                    %                             end
                end
                %end
            end
            
            legend('coh_{prev}<0','coh_{prev}>0')
            xlim([0,1])
            ylabel('Change in belief')
            %title('prev trial coh')
            %ht(1) = title('Data');
            
            xlabel('Current belief')
            
            p.next();
            
            SS = repmat(S,nsims,1);
            nn = length(repcoh(:));
            repcohprev = repmat(coh_prev,1,nsims);
            
            %cats = [repcohprev(:)<0,repcohprev(:)==0,repcohprev(:)>0];
            cats = [repcohprev(:)<0, repcohprev(:)>0];
            %cats = [repcohprev(:)==-0.256,repcohprev(:)==0,repcohprev(:)==0.256];
            
            [~,~,~,v] = index_prctile(belief_model_shifted(:)+randn(nn,1)*0.00001,[0:20:100]);
            for i=1:size(cats,2)
                %                     I = repcoh(:)==0 & cats(:,i)==1 & ~isnan(deltab_model(:)) & SS & ...
                %                          repmat(trnum_eff,nsims,1)==2;
                I = repcoh(:)==0 & cats(:,i)==1 & ~isnan(deltab_model(:)) & SS;
                
                [tt,xx,ss] = curva_media(deltab_model(:),v,I,0);
                terrorbar(tt,xx,ss,'marker','.','linestyle','-','color',colores(i,:),'linewidth',lw);
                hold all
                %                     [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),belief_model_shifted(:),npoints * nsims,0,I);
                %                     hold all
                %                     plot(tt{1},xx{1},'color',colores(i,:),'linestyle',lsty,'linewidth',lw);
            end
            
            
            xlabel('Current belief')
            
            
            same_ylim(p.h_ax);
            p.format('FontSize',14);
            
        end
        
        
        
        function combined_supp_fig3(obj)
            
            use_analog_conf = 0;
            
            
            p = publish_plot(3,4);
            set(gcf,'Position',[554  556  873  592])
            
            pacc = obj.acc_unsigned();
            for i=1:length(pacc.h_ax)
                pacc.current_ax(i);
                xlabel('');
            end
            if use_analog_conf
                pconf = obj.conf_continuous_unsigned_per_suj();
            else
                pconf = obj.conf_unsigned_per_suj();
            end
            for i=1:4
                p.copy_from_ax(pacc.h_ax(i),i);
                p.copy_from_ax(pconf.h_ax(i),i+4);
                p.copy_from_ax(pconf.h_ax(i+4),i+8);
            end
            
            p.unlabel_center_plots();
            close(pacc.h_fig);
            close(pconf.h_fig);
            BreakXAxis(p.h_ax,0.18,'break_scaling',4);
            set(p.h_ax,'xminortick','off','yminortick','off');
            
        end
        
        function p = combined_fig6_alt2(obj)
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            correct_only = 1;
            
            sbelief = belief;
            deltab = deltab_full;
            sbelief(choice==0) = 1-sbelief(choice==0);
            deltab(choice==0) = -1*deltab(choice==0);
            if correct_only
                deltab(correct==0) = nan;
            end
            
            deltab_model = deltab_model_full;
            if correct_only
                deltab_model(correct_model==0) = nan; % only correct trials
            end
            
            sbelief_model = belief_model;
            sbelief_model(choice_model==0) = 1-sbelief_model(choice_model==0);
            deltab_model(choice_model==0)  = -1*deltab_model(choice_model==0);
            
            %deltab_model = nanmean(deltab_model,2);
            
            p = publish_plot(2,2);
            set(gcf,'Position',[443  356  557  442])
            p.displace_ax([1,2],-0.08,2);
            
            npoints = 300;
            
            p.current_ax(1);
            % plot
            
            ig = 1;
            S = ismember(group,obj.G{ig});
            
            for i=1:length(uacoh)
                I = acoh==uacoh(i) & ~isnan(deltab) & S;
                [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab,ones(size(coh)),sbelief,npoints,0,I);
                hold all
                plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:))
            end
            
            ylabel({'DATA','\Deltabelief'})
            
            
            % model
            p.current_ax(3);
            
            SS = repmat(S,nsims,1);
            nn = length(repcoh(:));
            for i=1:length(uacoh)
                
                I = abs(repcoh(:))==uacoh(i) & ~isnan(deltab_model(:)) & SS;
                [tt,xx] = dotsanalysis.choice_against_duration_moving(deltab_model(:),ones(nn,1),sbelief_model(:),npoints * nsims,0,I);
                hold all
                plot(tt{1},xx{1},'color',obj.colores_uscoh(i,:))
            end
            
            xlabel({'Belief','(relative to choice)'})
            
            % now the split ones
            belief_rel_choice = belief;
            belief_rel_choice(choice==0) = 1-belief(choice==0);
            
            belief_model_rel_choice = belief_model;
            belief_model_rel_choice(choice_model==0) = 1-belief_model(choice_model==0);
            
            % data
            p.current_ax(2)
            I = belief_rel_choice<0.5 & S & ~isnan(deltab);
            dotsanalysis.plot_log(acoh(I),deltab(I),acoh(I),deltab(I),'color','r','break_scaling',0,'ylim',[0,0.08],'marker','o',...
                'markersize',8,'markerfacecolor','r','markeredgecolor','w');
            
            I = belief_rel_choice>0.5 & S & ~isnan(deltab);
            dotsanalysis.plot_log(acoh(I),deltab(I),acoh(I),deltab(I),'color','b','break_scaling',0,'ylim',[0,0.08],'marker','o',...
                'markersize',8,'markerfacecolor','b','markeredgecolor','w');
            
            
            p.current_ax(4);
            I = belief_model_rel_choice(:)<0.5 & SS(:) & ~isnan(deltab_model(:));
            hplot(1) = dotsanalysis.plot_log(abs(repcoh(I)),deltab_model(I),abs(repcoh(I)),deltab_model(I),'color','r','break_scaling',0,'ylim',[0,0.08],'marker','o',...
                'markersize',8,'markerfacecolor','r','markeredgecolor','w');
            
            I = belief_model_rel_choice(:)>0.5 & SS(:) & ~isnan(deltab_model(:));
            hplot(2) = dotsanalysis.plot_log(abs(repcoh(I)),deltab_model(I),abs(repcoh(I)),deltab_model(I),'color','b','break_scaling',0,'ylim',[0,0.08],'marker','o',...
                'markersize',8,'markerfacecolor','b','markeredgecolor','w');
            
            p.current_ax(3);
            ylabel('\Delta belief, MODEL')
            hl = legend(p.h_ax(4),hplot,'<0.5','>0.5');
            set(hl,'box','off');
            set(get(hl,'title'),'String',{'Belief','rel. choice'},'fontweight','normal')
            set(hl,'location','best');
            
            p.current_ax(4);
            %ylabel('\Delta belief, MODEL')
            xlabel({'Motion strength (%coh)',' '})
            
            set(p.h_ax(1:2),'xticklabel','');
            same_ylim(p.h_ax(1:2));
            same_ylim(p.h_ax(3:4));
            set(p.h_ax([2,4]),'ycolor','none');
            p.format('FontSize',14)
            set(hl,'FontSize',12);
            
            p.current_ax(1);
            plot([.5,.5],ylim,'k--')
            p.current_ax(3);
            plot([.5,.5],ylim,'k--')
            
            % BreakXAxis(p.h_ax([2,4]),0.022,'break_scaling',3);
            BreakXAxis(p.h_ax([2,4]),0.2,'break_scaling',7);
            set(p.h_ax,'xminortick','off');
            
            
            
        end
        
        
        
        function p = confidence_correct_vs_duration_split_by_coh(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            npoints = 15;
            
            p = publish_plot(1,nSujs+1);
            for ig=1:(nSujs+1)
                
                p.next();
                %                 I = ismember(repgroup,obj.G{ig}) & correct_model==1;
                %                 conf_aux = double(conf_high_model);
                %                 conf_aux(I==0) = nan;
                %                 [x,Y] = dotsanalysis.choice_against_duration(conf_aux(:),repcoh(:),repdotdur(:),npoints,0,I(:));
                
                
                % smooth
                if exist('smooth_model_high_given_correct_unsigned_d','var')
                    I = ismember(group,obj.G{ig});
                    mm = squeeze(nanmean(smooth_model_high_given_correct_unsigned_d(I,:,:)));
                    for i=1:size(mm,1)
                        plot(udurs,mm(i,:),'-','color',obj.colores_uscoh(i,:),'LineWidth',1)
                        hold all
                    end
                else
                    I = ismember(repgroup,obj.G{ig}) & correct_model==1;
                    conf_aux = double(conf_high_model);
                    conf_aux(I==0) = nan;
                    [x,Y] = dotsanalysis.choice_against_duration(nanmean(conf_aux,2),abs(coh),dotdur,npoints,0);
                    for i=1:size(Y,2)
                        plot(x,Y(:,i),'-','color',obj.colores_uscoh(i,:),'LineWidth',1)
                        hold all
                    end
                end
                
                
                I = ismember(group,obj.G{ig}) & correct==1;
                [x,Y,S] = dotsanalysis.choice_against_duration(high_conf,coh,dotdur,npoints,0,I);
                for i=1:size(Y,2)
                    terrorbar(x,Y(:,i),S(:,i),'marker','o','color',obj.colores_uscoh(i,:),'markerfacecolor',obj.colores_uscoh(i,:),...
                        'markeredgecolor','w','markersize',7,'LineStyle','--')
                    hold all
                end
                
                title(obj.titles{ig},'interpreter','none')
            end
            
            p.current_ax(1);
            xlabel('Stimulus duration')
            ylabel({'p high confidence',', correct trials'})
            
            set(p.h_ax,'xlim',[0.05,0.9],'tickdir','out')
            set(p.h_ax(2:end),'yticklabel','');
            p.format();
            %set(gcf,'Position',[43   240  1829   258])
            set(gcf,'Position',obj.figpos.onen);
            obj.figures_handle.confidence_correct_vs_duration = gcf;
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_conf_vs_dur'));
            %pconf = p;
            
        end
        
        
        function p = accuracy_vs_duration(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            p = publish_plot(1,nSujs+1);
            for ig=1:(nSujs+1)
                I = ismember(group,obj.G{ig});
                % I = ismember(group,G{ig}) & certain_its_biased;
                p.next();
                npoints = 15;
                
                if exist('smooth_model_correct_unsigned_d','var')
                    % smooth model
                    mm = squeeze(nanmean(smooth_model_correct_unsigned_d(I,:,:)));
                    for i=1:size(mm,1)
                        plot(udurs,mm(i,:),'-','color',obj.colores_uscoh(i,:),'LineWidth',1)
                        hold all
                    end
                else
                    [x,Y] = dotsanalysis.choice_against_duration(nanmean(correct_model,2),coh,dotdur,npoints,0,I);
                    for i=1:size(Y,2)
                        plot(x,Y(:,i),'-','color',obj.colores_uscoh(i,:),'LineWidth',1)
                        hold all
                    end
                end
                
                
                
                if 1
                    [x,Y,S] = dotsanalysis.choice_against_duration(correct,coh,dotdur,npoints,0,I);
                    for i=1:size(Y,2)
                        %plot(x,Y(:,i),'.-','color',obj.colores_uscoh(i,:))
                        terrorbar(x,Y(:,i),S(:,i),'marker','o','color',obj.colores_uscoh(i,:),'markerfacecolor',obj.colores_uscoh(i,:),...
                            'markeredgecolor','w','markersize',7,'LineStyle','--')
                        hold all
                    end
                end
                title(obj.titles{ig},'interpreter','none')
            end
            
            p.current_ax(1);
            xlabel('Stimulus duration (s)')
            ylabel('Accuracy')
            
            set(p.h_ax,'xlim',[0.05,0.9],'tickdir','out')
            set(p.h_ax(2:end),'yticklabel','');
            p.format();
            %set(gcf,'Position',[43   240  1829   258])
            set(gcf,'Position',obj.figpos.onen)
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_acc_vs_dur'));
            
            %pacc = p;
            obj.figures_handle.accuracy_vs_duration = gcf;
            
        end
        
        function p = belief_vs_trnum2(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            for i=1:length(uni_block)
                I = block_prior==uni_block(i) & trnum<=42;
                [xx_data{i},yy_data{i},ss_data{i}] = curva_media(belief,trnum,I,0);
                [xx_model{i},yy_model{i}] = curva_media_repmat(belief_model,trnum,I,0);
            end
            
            %[tt,xx,ss] = curva_media(belief ,trnum_eff,I,0);
            
            p = publish_plot(1,1);
            for i=1:length(uni_block)
                [~,hnice(i)] = niceBars2(xx_data{i}, yy_data{i}, ss_data{i}, obj.colores_us(i,:), 1);%.7
                hold all
            end
            set(hnice,'linestyle','none');
            for i=1:length(uni_block)
                plot(xx_model{i},yy_model{i},'color',obj.colores_us(i,:),'LineWidth',1);
                %                 terrorbar(xx_data{i},yy_data{i},ss_data{i},'color',obj.colores_us(i,:),'marker','o','markersize',5,...
                %                     'markerfacecolor',obj.colores_us(i,:),'markeredgecolor','w','linestyle','--');
            end
            xlabel('Trial number')
            ylabel('Belief')
            xlim([0,42]);
            
        end
        
        function belief_scatter_model_vs_data(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            %% do  Ss get the belief right?
            % get last belief
            
            p = publish_plot(1,nSujs+1);
            for i=1:length(obj.G)
                p.next();
                I = ismember(group,obj.G{i}) & trnum_back_eff==1;
                bm = nanmean(belief_model(I,:),2);
                plot(bm,belief(I),'linestyle','none','marker','o','markersize',6,'markerfacecolor','w','markeredgecolor','k')
                hold on
                plot([0.5,0.5],[0,1],'k--');
                plot([0,1],[0.5,0.5],'k--');
                title(obj.titles{i},'interpreter','none');
                if i==1
                    ylabel('End-belief, DATA')
                end
                xlabel('End-belief, MODEL')
                
                rsquared = RSquared(belief(I),bm);
                str = roundstr(rsquared,2);
                text(0.7,0.1,['R^2: ',str],'color',0.3*[1,1,1]);
            end
            
            set(p.h_ax(2:end),'yticklabel','');
            p.format();
            set(gcf,'Position',obj.figpos.onen)
            obj.figures_handle.end_belief = gcf;
            
        end
        
        function p = combined_fig6(obj,flag_average_simulations)
            
            
            % FIG for slide
            %             p = publish_plot(4,4);
            %             set(gcf,'Position',[966  853  594  485])
            %             rng(5861,'twister')
            %             % rng(5868,'twister')
            %             ses = randsample(uni_session,16);
            %             rng('shuffle','twister')
            %             for i=1:length(ses)
            %                 p.next;
            %                 I = session==ses(i);
            %                 plot(trnum_eff(I),belief_model(I,1:10),'color',0.7*[1,1,1],'LineWidth',0.4)
            %                 hold all
            %
            %                 plot(trnum_eff(I),belief(I),'k','LineWidth',2)
            %
            %                 area([0,42],[1,1],'facecolor','none','edgecolor','k','linewidth',0.332)
            %                 %plot([0,40],[.5,.5],'k--')
            %                 xlabel('Trial number')
            %                 ylabel('Belief')
            %             end
            %
            %             set(p.h_ax,'xlim',[0,42],'xtick',[5,20,35]);
            %
            %             p.unlabel_center_plots()
            %
            %             xlabel('Trial number')
            %             ylabel('Belief')
            
            
            plottype_flag = 2;
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            p = publish_plot(6,3);
            set(gcf,'Position',[631  171  353  627])
            p.combine([2,3,5,6]);
            p.combine([2,3,5,6]+6);
            p.combine([2,3,5,6]+12);
            p.clean_hax();
            p.shrink([7,8,9],0.8,[],1);
            p.displace_ax([1:6],-0.06,1)
            p.shrink([1:6],0.9,1)
            p.resize_vertical(1:6,ones(1,6),0.025);
            
            %p.number_the_plots();
            
            %duplicate the 9th
            p.new_axes('Position',get(p.h_ax(9),'Position'));
            %linkprop(p.h_ax([9,10]),'xlim','ylim');
            %linkprop(p.h_ax([1:6]),'xlim','ylim');
            
            
            %rng(58829,'twister');
            rng(58829+8,'twister');
            ses = randsample(uni_session,6);
            rng('shuffle','twister')
            colores_h = rainbow_colors(6,'colorType',6);
            for i=1:6
                p.current_ax(i);
                I = session==ses(i);
                plot(trnum_eff(I),belief_model(I,1:10),'color',0.7*[1,1,1],'LineWidth',0.4);
                hold all
                
                hp(i) = plot(trnum_eff(I),belief(I),'k','LineWidth',2);
                
                ha(i) = area([0,42],[1,1],'facecolor','none','edgecolor','k','linewidth',0.332);
                hold on
                
                switch plottype_flag
                    case 1
                        plot(35,0.8,'marker','o','markersize',6,'markerfacecolor',colores_h(i,:),'markeredgecolor',colores_h(i,:));
                    case 2
                        set(hp(i),'color',colores_h(i,:));
                    case 3
                        set(gca,'ycolor',colores_h(i,:));
                    case 4
                        set(ha(i),'edgecolor',colores_h(i,:));
                        
                end
            end
            
            xlabel('Trial number')
            ylabel('Belief')
            
            
            p.current_ax(7); %belief vs block bias, model and data
            for i=1:length(uni_block)
                I = block_prior==uni_block(i) & trnum<=42;
                [xx_data{i},yy_data{i},ss_data{i}] = curva_media(belief,trnum,I,0);
                [xx_model{i},yy_model{i}] = curva_media_repmat(belief_model,trnum,I,0);
            end
            
            %[tt,xx,ss] = curva_media(belief ,trnum_eff,I,0);
            
            
            for i=1:length(uni_block)
                [~,hnice(i)] = niceBars2(xx_data{i}, yy_data{i}, ss_data{i}, obj.colores_us(i,:), 1);%.7
                hold all
            end
            set(hnice,'linestyle','none');
            for i=1:length(uni_block)
                plot(xx_model{i},yy_model{i},'color',obj.colores_us(i,:),'LineWidth',1);
                %                 terrorbar(xx_data{i},yy_data{i},ss_data{i},'color',obj.colores_us(i,:),'marker','o','markersize',5,...
                %                     'markerfacecolor',obj.colores_us(i,:),'markeredgecolor','w','linestyle','--');
            end
            xlabel('Trial number')
            ylabel('Belief')
            xlim([0,42]);
            
            
            p.current_ax(8); % end belief, model and data
            
            %             area([0,0.5],[0,0],.5,'facecolor',0.8*[1,1,1],'edgecolor','none')
            %             hold all
            %             area([0.5,1],[1,1],.5,'facecolor',0.8*[1,1,1],'edgecolor','none')
            
            I = ismember(group,obj.G{1}) & trnum_back_eff==1;
            
            if flag_average_simulations
                bm = nanmean(belief_model,2); %version in submitted paper
                rsquared = RSquared(belief(I),bm(I));
            else
                bm = nanmean(belief_model(:,1),2);
                for k=1:size(belief_model,2)
                    v_rsquared(k) = RSquared(belief(I),belief_model(I,k));
                end
                rsquared = mean(v_rsquared);
            end
            plot(bm(I),belief(I),'linestyle','none','marker','o','markersize',6,'markerfacecolor','k','markeredgecolor','w','LineWidth',0.33)
            
            hold on
            %I = find(ismember(group,obj.G{1}) & trnum_back_eff==1 & ismember(session,ses));
            for j=1:length(ses)
                I = session==ses(j) & trnum_back_eff == 1;
                plot(bm(I),belief(I),'linestyle','none','marker','o','markersize',7,...
                    'markerfacecolor',colores_h(j,:),'markeredgecolor',colores_h(j,:))
            end
            %plot([0.5,0.5],[0,1],'k--');
            %plot([0,1],[0.5,0.5],'k--');
            xlabel('End-of-block belief, Prediction')
            ylabel('End-of-block belief, Data')
            
            str = roundstr(rsquared,2);
            if flag_average_simulations
                text(0.8,0.1,['R^2: ',str],'color',0.3*[1,1,1]);
            else
                text(0.8,0.1,['\langleR^2\rangle: ',str],'color',0.3*[1,1,1]);
            end
            %axis square
            
            p.current_ax(9); % prop correct and conf correct
            
            isLastTrial = trnum_back_eff==1;
            
            %[tt,xx,ss] = curva_media(belief_conf,block_prior_rel,belief_correct==0 & I,0);
            %terrorbar(tt,xx,ss,'color','r','marker','o','markersize',7,'markerfacecolor','r','markeredgecolor','w')
            
            
            
            option_flag = 3;
            
            switch option_flag
                case 1
                    
                    belief_conf = max(belief,1-belief);
                    belief_conf_model = max(belief_model,1-belief_model);
                    
                    belief_correct  = belief>0.5;
                    belief_correct(block_prior<0.5) = 1-belief_correct(block_prior<0.5);
                    
                    belief_correct_model = belief_model>0.5;
                    I = block_prior<0.5;
                    belief_correct_model(I,:) = 1-belief_correct_model(I,:);
                    
                    I = ismember(group,obj.G{1}) & isLastTrial==1;
                    Ir = repmat(I,1,nsims);
                    
                    [tt,xx,ss] = curva_media_repmat(belief_correct_model,block_prior_rel,I,0);
                    h(1) = plot(tt,xx,'color','k','linestyle','-');
                    hold all
                    [tt,xx,ss] = curva_media(belief_correct,block_prior_rel,I,0);
                    h(2) = terrorbar(tt,xx,ss,'color','k','marker','o','markersize',7,'markerfacecolor','k','markeredgecolor','w','linestyle','none');
                    ylabel({'P correct about','sign of bias'})
                    
                    hl = legend(h(2:-1:1),'data','prediction');
                    set(hl,'location','best','box','on');
                    
                    
                    p.current_ax(10); % prop correct and conf correct
                    set(gca,'color','none','yaxislocation','right','ycolor','b');
                    
                    hold all
                    
                    bpr = repmat(block_prior_rel,1,nsims);
                    
                    [tt,xx,ss] = curva_media(belief_conf_model(:),bpr(:),belief_correct_model(:) & Ir(:),0);
                    plot(tt,xx,'color','r','linestyle','-');
                    
                    [tt,xx,ss] = curva_media(belief_conf,block_prior_rel,belief_correct & I,0);
                    terrorbar(tt,xx,ss,'color','r','marker','o','markersize',7,'markerfacecolor','r','markeredgecolor','w','linestyle','none')
                    
                    ylabel('Belief strength');
                    xlabel('Bias strength');
                    
                case 2
                    
                    set(p.h_ax(10),'visible','off');
                    p.current_ax(9);
                    
                    I = ismember(group,obj.G{1}) & isLastTrial==1;
                    Ir = repmat(I,1,nsims);
                    
                    bpr = repmat(block_prior,1,nsims);
                    
                    [tt,xx,ss] = curva_media(belief_model(:),bpr(:),Ir(:),0);
                    plot(tt,xx,'color','k','linestyle','-');
                    
                    hold all
                    
                    [tt,xx,ss] = curva_media(belief, block_prior, I, 0);
                    terrorbar(tt,xx,ss,'color','k','marker','o','markersize',7,'markerfacecolor','k','markeredgecolor','w','linestyle','none')
                    
                    ylabel('Belief');
                    %ylabel('End-of-block belief')
                    xlabel('Base-rate of block');
                    
                    
                    set(gca,'xtick',0:0.2:1);
                    
                case 3
                    
                    set(p.h_ax(10),'visible','off');
                    p.current_ax(9);
                    p.shrink(9,0.9,0.9);
                    
                    I = ismember(group,obj.G{1}) & isLastTrial==1;
                    
                    Ir = repmat(I,1,nsims);
                    bpr = repmat(block_prior,1,nsims);
                    [ttm,xxm,ssm] = curva_media(belief_model(:),bpr(:),Ir(:),0);
                    
                    %                     if flag_average_simulations
                    %                         Ir = repmat(I,1,nsims);
                    %                         bpr = repmat(block_prior,1,nsims);
                    %                         [ttm,xxm,ssm] = curva_media(belief_model(:),bpr(:),Ir(:),0);
                    %                     else % only one simulation
                    %                         [ttm,xxm,ssm] = curva_media(belief_model(:,1),block_prior(:,1),I,0);
                    %                     end
                    
                    %plot(tt,xx,'color','k','linestyle','-');
                    
                    hold all
                    
                    [ttd,xxd,ssd] = curva_media(belief, block_prior, I, 0);
                    %terrorbar(tt,xx,ss,'color','k','marker','o','markersize',7,'markerfacecolor','k','markeredgecolor','w','linestyle','none')
                    
                    
                    plot([0,1],[0,1],'k--')
                    hold all
                    
                    for i=1:length(xxm)
                        colo = obj.colores_us(i,:);
                        plot(xxm(i),xxd(i),'o','markersize',9,'markerfacecolor',colo,'markeredgecolor',colo,'linestyle','none');
                        plot([xxm(i)-ssm(i),xxm(i)+ssm(i)],[xxd(i),xxd(i)],'color',colo);
                        plot([xxm(i),xxm(i)],[xxd(i)-ssd(i),xxd(i)+ssd(i)],'color',colo);
                    end
                    
                    xlabel('End-of-block belief, Prediction')
                    ylabel('End-of-block belief, Data')
                    
                    %set(gca,'xtick',0:0.2:1);
                    
            end
            %[tt,xx,ss] = curva_media(belief_conf_model(:),bpr(:),belief_correct_model(:)==0 & Ir(:),0);
            %terrorbar(tt,xx,ss,'color','r','linestyle','-');
            
            %hl = legend('correct, data','error, data','correct, model','error, model');
            %set(hl,'location','best','fontsize',10,'box','on');
            
            
            % ylabel({'Confidence in ','bias judgement'})
            
            
            set(p.h_ax(1:6),'xlim',[0,42],'xtick',[5,20,35]);
            set(p.h_ax(1:6),'ylim',[0,1])
            set(p.h_ax(1:5),'xticklabel','','xcolor','none');
            
            if option_flag==1
                set(p.h_ax(9:10),'ylim',[0.6,1]);
                set(p.h_ax(9),'color','none');
                set(p.h_ax(10),'yaxislocation','right','color','none','ycolor','r')
            else
                
                
            end
            
            set(p.h_ax,'tickdir','out');
            p.format('FontSize',11);
            
            p.current_ax(9);
            nontouching_spines(gca);
            %             nontouching_spines(p.h_ax(1:6));
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% UNSIGNED PLOT %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % keep
        function p = acc_unsigned(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            % u = unique(dat(1).block_prior);
            
            u = unique(block_bias);
            
            cmodela = obj.m.smooth_model_correct_unsigned;
            
            p = publish_plot(1,nSujs+1);
            for ig=1:(nSujs+1)
                I = ismember(group,obj.G{ig});
                p.next();
                for i=1:length(u)
                    J = block_bias==u(i) & I;
                    ntr(i) = sum(J);
                    dotsanalysis.plot_log(coh(J),correct(J),[],[],'linestyle','none','color',color_block_bias(i,:),'marker','o','markersize',6,'break_scaling',0)
                    hold all
                    if smooth_flag
                        
                        %             JJ = block_bias==u(i) & I;
                        dotsanalysis.plot_log(aucohs,nanmean(cmodela(J,:)),[],[],'linestyle','-','color',color_block_bias(i,:),'marker','none','markerfacecolor','w',...
                            'markersize',6,'break_scaling',0,'show_errors',0);
                    else
                        cmodel = nanmean(pcorrect_model,2);
                        dotsanalysis.plot_log(coh(J),cmodel(J),[],[],'linestyle','-','color',color_block_bias(i,:),'marker','none','markerfacecolor','w',...
                            'markersize',6,'break_scaling',0,'show_errors',0)
                    end
                end
                %text(0.4,0.5,['N=',num2str(sum(ntr))]);
                if ig==1
                    %         legend_n(u,'hline',h)
                    xlabel('Coherence')
                    ylabel('p correct')
                else
                    xlabel(' ') % I hate matlab
                end
                title(obj.titles{ig},'interpreter','none')
            end
            
            % set(p.h_ax,'xlim',[-0.6,0.6])
            set(p.h_ax(2:end),'yticklabel','');
            p.format('FontSize',14);
            
            set(p.h_ax,'tickdir','out')
            
            set(gcf,'Position',obj.figpos.onen_u);
            %set(gcf,'Position',[1   578  1440   184])
            
            
            pause(1);
            % BreakXAxis(p.h_ax,0.022,'break_scaling',10);
            BreakXAxis(p.h_ax,0.18,'break_scaling',10);
            
            set(p.h_ax,'XMinorTick','off','YMinorTick','off');
            
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_acc_unsigned'));
            obj.figures_handle.acc_unsigned = gcf;
            %phandle_pacc = p;
            
            
            
            % 2nd fig: fig_acc_unsigned_per_stage
            
            
            %% per stage
            % tr = {1:5,6:10,11:15};
            % tr_split = {1:10,11:20,21:30};
            
            %
        end
        
        
        % keep
        function p = conf_unsigned_per_suj(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            %% confidence
            % bl = round(100*(block_prior.*choice+(1-block_prior).*(1-choice)))/100;%prior rel to choice
            
            
            % aux calc min trials:
            conditions = [abs(coh),group,bl];
            [~,~,uni_conditions,~,idx_cond] = average_per_condition(ones(size(coh)), conditions,'filter',correct==0);
            
            a = [Rtable(idx_cond(~isnan(idx_cond))),uni_conditions];
            
            
            
            min_num_trials_errors = 10;
            first_coh_to_exclude = a(find(a(:,1)<min_num_trials_errors,1),2); % all sujs together
            
            for i=1:nSujs
                ss = obj.uni_sujs(i);
                I = a(:,3)==ss;
                first_coh_to_exclude(i+1) = a(find(a(:,1)<min_num_trials_errors & I,1),2); %per suj
            end
            
            %
            
            p = publish_plot(2,nSujs+1);
            % set(gcf,'Position',[586   28  307  777])
            %set(gcf,'Position',[74   452  1881   441])
            set(gcf,'Position',obj.figpos.twon_u);
            
            for ig=1:(nSujs+1)
                F = ismember(group,obj.G{ig});
                p.current_ax(ig);
                
                clear ntr
                for i=1:length(uni_prior_rel_choice)
                    
                    I = bl==uni_prior_rel_choice(i) & F;
                    
                    debug_flag = 0;
                    
                    if debug_flag
                        J = I & correct_model(:,1)==1;
                        xdata = abs(coh(J));
                        ydata = conf_high_model(J,1)==1;
                        
                    else
                        J = I & correct==1;
                        xdata = abs(coh(J));
                        ydata = high_conf(J)==1;
                        ntr(i) = nansum(J);
                    end
                    
                    
                    if smooth_flag
                        xmodel = aucohs;
                        ymodel = nanmean(obj.m.smooth_model_high_given_correct_unsigned(I,:));
                        %                         ymodel = nanmean(p_high_given_correct(J,:));
                        %                         xmodel = ucohs_zplus;
                    else
                        conf_model_aux = double(conf_high_model);
                        conf_model_aux(correct_model==0) = nan;
                        [xmodel,ymodel] = curva_media_repmat(conf_model_aux,abs(coh),I);
                    end
                    dotsanalysis.plot_log(xdata,ydata,xmodel,ymodel,'color',coloresm(i,:),'LineStyle','-','ylim',[0.05,1],...
                        'marker','o','markersize',5,'break_scaling',0)
                    %         ylabel('P high confidence');
                    %         xlabel('Motion strength (%coh)')
                end
                
                %ht(ig,1) = text(0.4,0.2,['N=',num2str(sum(ntr))]);
                
                p.current_ax(ig+nSujs+1);
                
                clear ntr
                for i=1:length(uni_prior_rel_choice)
                    
                    I = bl==uni_prior_rel_choice(i) & F;
                    
                    J = I & correct==0 & abs(coh)<first_coh_to_exclude(ig);
                    
                    if debug_flag
                        J = I & correct_model(:,1)==0 & abs(coh)<first_coh_to_exclude(ig);
                        xdata = abs(coh(J));
                        ydata = conf_high_model(J)==1;
                        
                    else
                        J = I & correct==0 & abs(coh)<first_coh_to_exclude(ig);
                        xdata = abs(coh(J));
                        ydata = high_conf(J)==1;
                        ntr(i) = nansum(J);
                    end
                    
                    last_coh_to_include = max(uacoh(uacoh<first_coh_to_exclude(ig)));
                    
                    ff = aucohs<=(last_coh_to_include+0.02);
                    if smooth_flag
                        xmodel = aucohs(ff);
                        ymodel = nanmean(obj.m.smooth_model_high_given_error_unsigned(I,ff));
                        %                         ymodel = nanmean(p_high_given_error(J,ff));
                        %                         xmodel = ucohs_zplus(ff);
                    else
                        conf_model_aux = double(conf_high_model);
                        conf_model_aux(correct_model==1) = nan;
                        [xmodel,ymodel] = curva_media_repmat(conf_model_aux,abs(coh),I);
                    end
                    dotsanalysis.plot_log(xdata,ydata,xmodel,ymodel,'color',coloresm(i,:),'LineStyle','-','ylim',[0.05,1],...
                        'marker','o','markersize',5,'break_scaling',0)
                    %         ylabel('P high confidence');
                    %         xlabel('Motion strength (%coh)')
                end
                
                %ht(ig,1) = text(0.4,0.2,['N=',num2str(sum(ntr))]);
                
            end
            
            p.current_ax(1);
            %title('Correct')
            ylabel({'P high confidence','correct trials'})
            
            p.current_ax(nSujs+2);
            %title('Errors')
            ylabel({'P high confidence','error trials'})
            
            p.current_ax(nSujs+2);
            xlabel('Motion strength (%coh)')
            % ylabel({'P high confidence','error trials'});
            
            p.current_ax(1);
            % ylabel({'P high confidence','correct trials'});
            
            set(p.h_ax,'XMinorTick','off','YMinorTick','off');
            p.format('FontSize',14);
            set(p.h_ax,'ytick',0.2:0.2:1,'tickdir','out')
            set(p.h_ax(2:1:nSujs+1),'yticklabel','');
            set(p.h_ax((3+nSujs):end),'yticklabel','');
            %set(p.h_ax(1:nSujs+1),'xticklabel','');
            
            pause(.1);
            % BreakXAxis(p.h_ax,0.022,'break_scaling',5);
            BreakXAxis(p.h_ax,0.18,'break_scaling',7);
            
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_acc_unsigned_conf_per_acc_per_suj'));
            obj.figures_handle.acc_unsigned_conf_per_acc_per_suj_2 = gcf;
            
            drawnow
            
        end
        
        % keep
        function p = conf_continuous_unsigned_per_suj(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            %% confidence
            % bl = round(100*(block_prior.*choice+(1-block_prior).*(1-choice)))/100;%prior rel to choice
            
            
            
            if exist('conf_remapped','var')
                confi = conf_remapped;
            else
                confi = conf/2+0.5;
            end
            
            % aux calc min trials:
            conditions = [abs(coh),group,bl];
            [~,~,uni_conditions,~,idx_cond] = average_per_condition(ones(size(coh)), conditions,'filter',correct==0);
            
            a = [Rtable(idx_cond(~isnan(idx_cond))),uni_conditions];
            
            
            
            min_num_trials_errors = 10;
            first_coh_to_exclude = a(find(a(:,1)<min_num_trials_errors,1),2); % all sujs together
            
            for i=1:nSujs
                ss = obj.uni_sujs(i);
                I = a(:,3)==ss;
                first_coh_to_exclude(i+1) = a(find(a(:,1)<min_num_trials_errors & I,1),2); %per suj
            end
            
            %
            
            p = publish_plot(2,nSujs+1);
            % set(gcf,'Position',[586   28  307  777])
            %set(gcf,'Position',[74   452  1881   441])
            set(gcf,'Position',obj.figpos.twon_u);
            
            for ig=1:(nSujs+1)
                F = ismember(group,obj.G{ig});
                p.current_ax(ig);
                
                clear ntr
                for i=1:length(uni_prior_rel_choice)
                    
                    I = bl==uni_prior_rel_choice(i) & F;
                    
                    debug_flag = 0;
                    
                    if debug_flag
                        J = I & correct_model(:,1)==1;
                        xdata = abs(coh(J));
                        ydata = conf_model(J,1);
                        
                    else
                        J = I & correct==1;
                        xdata = abs(coh(J));
                        %ydata = high_conf(J)==1;
                        ydata = confi(J);
                    end
                    ntr(i) = nansum(J);
                    
                    if isfield(obj.m,'smooth_model_econf_given_correct_unsigned') % to do:write the smooth part
                        xmodel = aucohs;
                        ymodel = nanmean(obj.m.smooth_model_econf_given_correct_unsigned(I,:));
                        %                         ymodel = nanmean(p_high_given_correct(J,:));
                        %                         xmodel = ucohs_zplus;
                    else
                        conf_model_aux = conf_model;
                        conf_model_aux(correct_model==0) = nan;
                        [xmodel,ymodel] = curva_media_repmat(conf_model_aux,abs(coh),I);
                    end
                    dotsanalysis.plot_log(xdata,ydata,xmodel,ymodel,'color',coloresm(i,:),'LineStyle','-','ylim',[0.5,1],...
                        'marker','o','markersize',5,'break_scaling',0)
                    %         ylabel('P high confidence');
                    %         xlabel('Motion strength (%coh)')
                end
                
                %ht(ig,1) = text(0.3,0.6,['N=',num2str(sum(ntr))]);
                
                p.current_ax(ig+nSujs+1);
                
                clear ntr
                for i=1:length(uni_prior_rel_choice)
                    
                    I = bl==uni_prior_rel_choice(i) & F;
                    
                    J = I & correct==0 & abs(coh)<first_coh_to_exclude(ig);
                    
                    if debug_flag
                        J = I & correct_model(:,1)==0 & abs(coh)<first_coh_to_exclude(ig);
                        xdata = abs(coh(J));
                        %ydata = conf_high_model(J)==1;
                        ydata = conf_model(J,1);
                        
                    else
                        J = I & correct==0 & abs(coh)<first_coh_to_exclude(ig);
                        xdata = abs(coh(J));
                        ydata = confi(J);
                    end
                    ntr(i) = nansum(J);
                    
                    last_coh_to_include = max(uacoh(uacoh<first_coh_to_exclude(ig)));
                    
                    ff = aucohs<=(last_coh_to_include+0.02);
                    if isfield(obj.m,'smooth_model_econf_given_correct_unsigned') % to do:write the smooth part
                        xmodel = aucohs(ff);
                        ymodel = nanmean(obj.m.smooth_model_econf_given_error_unsigned(I,ff));
                        %ymodel = nanmean(obj.m.smooth_model_high_given_error_unsigned(I,ff));
                        %                         ymodel = nanmean(p_high_given_error(J,ff));
                        %                         xmodel = ucohs_zplus(ff);
                    else
                        conf_model_aux = conf_model;
                        conf_model_aux(correct_model==1) = nan;
                        [xmodel,ymodel] = curva_media_repmat(conf_model_aux,abs(coh),I);
                    end
                    dotsanalysis.plot_log(xdata,ydata,xmodel,ymodel,'color',coloresm(i,:),'LineStyle','-','ylim',[0.5,1],...
                        'marker','o','markersize',5,'break_scaling',0)
                    %         ylabel('P high confidence');
                    %         xlabel('Motion strength (%coh)')
                end
                
                %ht(ig,1) = text(0.3,0.6,['N=',num2str(sum(ntr))]);
                
            end
            
            set(p.h_ax,'color','none');
            p.current_ax(1);
            %title('Correct')
            ylabel({'Confidence','correct trials'})
            
            p.current_ax(nSujs+2);
            %title('Errors')
            ylabel({'Confidence','error trials'})
            
            p.current_ax(nSujs+2);
            xlabel('Motion strength (%coh)')
            % ylabel({'P high confidence','error trials'});
            
            p.current_ax(1);
            % ylabel({'P high confidence','correct trials'});
            
            set(p.h_ax,'XMinorTick','off','YMinorTick','off');
            p.format('FontSize',14);
            set(p.h_ax,'ytick',0.6:0.2:1,'tickdir','out')
            set(p.h_ax(2:1:nSujs+1),'yticklabel','');
            set(p.h_ax((3+nSujs):end),'yticklabel','');
            %set(p.h_ax(1:nSujs+1),'xticklabel','');
            
            pause(.1);
            % BreakXAxis(p.h_ax,0.022,'break_scaling',5);
            BreakXAxis(p.h_ax,0.18,'break_scaling',7);
            
            %export_fig(gcf,'-pdf',fullfile(folder,'figures','fig_acc_unsigned_conf_per_acc_per_suj'));
            obj.figures_handle.acc_unsigned_conf_continous_per_acc_per_suj = gcf;
            
            drawnow
            
        end
        
        % keep
        function plot_weights(obj,average_flag)
            
            if nargin==1
                average_flag=0;
            end
            
            %% plot weights, get values
            if isfield(obj.m,'K')
                
                p = publish_plot(1,1);
                switch average_flag
                    case 0
                        h = plot(obj.m.K',obj.m.prior_K_ini','.-');
                        legend_n(obj.uni_sujs,'title','subject');
                        
                    case 1
                        mK = nanmean(obj.m.K);
                        mPriorIni = nanmean(obj.m.prior_K_ini);
                        
                        bar(mK,mPriorIni);
                        p.new_axes('Position',get(p.h_ax,'Position'));
                        
                        linkaxes(p.h_ax);
                        
                        p.current_ax(1);
                        Kgen = 0:0.2:1;
                        PriorGen = ones(1,6)/6;
                        hb(1) = bar(Kgen,PriorGen);
                        
                        p.current_ax(2);
                        hb(2) = bar(mK,mPriorIni);
                        
                        same_ylim(p.h_ax);
                        same_xlim(p.h_ax);
                        set(p.h_ax,'color','none');
                        set(hb(1),'facecolor','w');
                        set(hb(2),'facecolor','k');
                        
                        
                        set(hb,'barwidth',0.7);
                        set(p.h_ax,'xtick',0:0.2:1,'tickdir','out','ticklength',[0.02,0.02],...
                            'ytick',[]);
                        
                        set(p.h_ax(1),'visible','off');
                        ylabel('$p_0(B)$','interpreter','latex');
                        set(gcf,'Position',[930  850  581  375])
                        
                    case 2
                        
                        p = publish_plot(1,obj.nsujs);
                        set(gcf,'Position',[164  356  946  223])
                        
                        for i=1:obj.nsujs
                            
                            K = obj.m.K(i,:);
                            pK = obj.m.prior_K_ini(i,:);
                            
                            p.next();
                            hb(i) = bar(K,pK);
                            
                            ht(i) = title(['Subject ',num2str(obj.uni_sujs(i))]);
                        end
                        
                        same_ylim(p.h_ax);
                        set(p.h_ax,'xlim',[-0.05,1.05]);
                        set(p.h_ax,'color','none');
                        
                        set(ht,'fontweight','normal');
                        
                        %set(hb,'barwidth',0.7);
                        set(p.h_ax,'xtick',0:0.2:1,'tickdir','out','ticklength',[0.02,0.02],...
                            'ytick',[]);
                        
                        ratio = diff(obj.m.K(:,1:2),[],2)./(min(diff(obj.m.K(:,1:2),[],2)));
                        for i=1:obj.nsujs
                            set(hb(i),'barwidth',0.8/ratio(i));
                        end
                        
                        set(hb,'facecolor','k');
                        
                        p.current_ax(1);
                        ylabel('$p_0(B)$','interpreter','latex');
                        
                        p.format('FontSize',15);
                        %set(gcf,'Position',[930  850  581  375])
                        
                    case 3
                        
                        p = publish_plot(1,obj.nsujs);
                        set(gcf,'Position',[164  356  946  223])
                        
                        for i=1:obj.nsujs
                            
                            K = obj.m.K(i,:);
                            pK = obj.m.prior_K_ini(i,:);
                            
                            p.next();
                            
                            hb(i) = stem(K,pK);
                            
                            ht(i) = title(['Subject ',num2str(obj.uni_sujs(i))]);
                        end
                        
                        same_ylim(p.h_ax);
                        set(p.h_ax,'xlim',[-0.05,1.05]);
                        set(p.h_ax,'color','none');
                        
                        set(ht,'fontweight','normal');
                        
                        %set(hb,'barwidth',0.7);
                        set(p.h_ax,'xtick',0:0.2:1,'tickdir','out','ticklength',[0.02,0.02],'ylim',[0,0.5]);
                        
                        %color_bars = [223,64,152]/255;
                        color_bars = [0,0,0];
                        set(hb,'color',color_bars,'marker','o','markeredgecolor',color_bars,'markerfacecolor',color_bars,'linestyle','-','markersize',10,...
                            'LineWidth',2);
                        
                        
                        p.current_ax(1);
                        ylabel('$p_0(B)$','interpreter','latex');
                        
                        p.format('FontSize',15);
                        
                        
                end
                
                xlabel('Block bias (B)')
                ylabel('$p_0(B)$','interpreter','latex');
                
                p.format();
                
                drawnow
                
            end
            
        end
        
        
        
        
        % keep
        function influence_bias_0coh_option2(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            
            %% influence of previous conf and choice on current: both matter?
            prev_belief = belief - deltab_full;
            
            cutoff = 0.95;
            I = ((prev_belief>cutoff) | (prev_belief<(1-cutoff))) & coh==0;
            choice_rel = choice.*(prev_belief>=0.5) + (1-choice).*(prev_belief<0.5);
            [ttd,xxd,ssd] = curva_media(choice_rel,block_prior_rel,I,0);
            
            
            prev_belief_model = belief_model - deltab_model_full;
            Im = ((prev_belief_model(:)>cutoff) | (prev_belief_model(:)<(1-cutoff))) & repcoh(:)==0;
            %rep_block_prior = repmat(block_prior,[1,nsims]);
            choice_model_rel = choice_model(:).*(prev_belief_model(:)>=0.5) + (1-choice_model(:)).*(prev_belief_model(:)<0.5);
            rep_block_prior_rel = repmat(block_prior_rel,nsims,1);
            [ttm,xxm,ssm] = curva_media(choice_model_rel(:),rep_block_prior_rel(:),Im,0);
            
            
            p = publish_plot(1,1);
            set(gcf,'Position',[936  760  531  430])
            
            p.shrink(1,0.8,.8,1,1);
            %             p.displace_ax(1,0.2,1);
            %             p.displace_ax(1,0.1,2);
            
            %colo = obj.colores_us(4:end,:);
            colo = movshon_colors(3);
            for i=1:length(xxm)
                plot(xxm(i),xxd(i),'marker','o','markerfacecolor',colo(i,:),'linestyle','none','markersize',11,'markeredgecolor',colo(i,:));
                h = errorbarxy(xxm(i),xxd(i),ssm(i), ssm(i), ssd(i), ssd(i));
                set(h,'color',colo(i,:),'linewidth',1);
                hold all
            end
            
            same_xylim(p.h_ax)
            hh = refline(1,0);
            %set(hh,'Color',[195,69,149]/255);
            set(hh,'Color','k','linewidth',1,'linestyle','--');
            
            xlabel('Prediction');
            ylabel('Data');
            
            %title({'p. choices in direction of bias','for 0% trials'});
            %hl = legend('Data','Model');
            %set(hl,'location','best','FontSize',13);
            
            %xlab = {{'0.4','or 0.6'},{'0.2','or 0.8'},{'0','or 1'}};
            %ht = cell_in_xticklabel(gca, 0.6:0.2:1, xlab);
            
            %set(p.h_ax,'xtick',[0.6:0.2:1],'xticklabel',xlab);
            
            p.format('FontSize',18);
            
            nontouching_spines(p.h_ax,'ticklength',0.01);
            
            %export_fig('-pdf',fullfile(folder,'figures','fig_influence_bias_0coh'));
            obj.figures_handle.influence_bias_0coh = gcf;
            
            
            % stats data
            depvar = choice_rel(I);
            prev_belief_rel = max(prev_belief,1-prev_belief);
            %indepvar = {'group',adummyvar(group(I)),'block_prior',block_prior_rel(I),'dotdur',dotdur(I),'prev_belief_model',prev_belief_rel(I)};
            indepvar = {'group',adummyvar(group(I)),'block_prior',block_prior_rel(I),'dotdur',dotdur(I)};
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
            stats.p(idx.block_prior)
            
            % stats model
            depvar = choice_model_rel(Im);
            %prev_belief_rel_model = max(prev_belief_model,1-prev_belief_model);
            %             indepvar = {'group',adummyvar(repgroup(Im)),'block_prior',rep_block_prior_rel(Im),'dotdur',repdotdur(Im),...
            %                 'prev_belief_model',prev_belief_rel_model(Im)};
            indepvar = {'group',adummyvar(repgroup(Im)),'block_prior',rep_block_prior_rel(Im),'dotdur',repdotdur(Im)};
            [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);
            stats.p(idx.block_prior)
            
            
            drawnow
            
        end
        
        % keep
        function p = confidence_transformation(obj)
            
            struct2vars(obj.d);
            struct2vars(obj.m);
            nSujs = obj.nsujs;
            
            conf_re = 0.5+0.5*conf;
            p = publish_plot(1,nSujs+1);
            for ig=1:(nSujs+1)
                I = ismember(group,obj.G{ig});
                p.next();
                plot(sort(conf_re(I)),sort(conf_model(I,1)),'.')
                
                
                % also plot accuracy in bins? moving average?
                acc_flag = 1;
                switch acc_flag
                    case 1
                        npoints = 150;
                        [X,Y] = dotsanalysis.choice_against_duration_moving(correct,ones(size(correct)),conf_re,npoints,0,I);
                        hold all;plot(X{1},Y{1},'color',[0,0.5,0],'LineWidth',2);
                    case 2
                        [~,~,~,v] = index_prctile(conf_re(I)+randn(sum(I),1)*0.00001,[0:5:100]);
                        [tt,xx,ss] = curva_media(correct(I),v,[],0);
                        hold all
                        plot(tt,xx,'marker','.','color',[0,0.5,0],'LineWidth',2);
                end
                
                %                 if nanmean(conf_re(I))>nanmean(conf_model(I,1))
                %                     text(0.6,0.45,'overconfident');
                %                 else
                %                     text(0.6,0.45,'underconfident');
                %                 end
                
                
                xlabel('Confidence, sorted');
                if ig==1
                    ylabel('Model conf, sorted');
                end
                title(obj.titles{ig},'interpreter','none');
                axis square
                h = refline(1,0);
                set(h,'LineWidth',1.5,'color','k','linestyle',':');
            end
            set(p.h_fig,'position',obj.figpos.onen);
            p.format();
            
        end
        
        
        drawnow
        
        
        
    end
end

% function hlines = calc_and_plot_kernels(signed_var,group,session,belief,trnum_eff,hax)
%
% [ev_coh,toff_coh,datamp,last_tr,last_belief] = calc_belief_kernels(signed_var,group,session,belief,trnum_eff);
%
% %I = block_prior(last_tr)~=0 & block_prior(last_tr)~=1;
%
% set(gcf,'CurrentAxes',hax(1));
% hlines(1) = plot(nanmean(ev_coh(last_belief>0.5,:)));
% hold all
% hlines(2) = plot(nanmean(ev_coh(last_belief<0.5,:)));
% ylabel('Average signed coherence')
% pxlim
%
% set(gcf,'CurrentAxes',hax(2));
% hlines(3) = plot(toff_coh,nanmean(datamp(:,last_belief>0.5),2));
% hold all
% hlines(4) = plot(toff_coh,nanmean(datamp(:,last_belief<0.5),2));
% pxlim
%
% end
