%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing data for TFA
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%steps:
%1) loading data with parameters specified in script "PARAMS"
%2) applying filters
%3) cleaning channel names
%4) checking triggers
%5) exploring signal channelwise, and exclude noisy channels
%6) applying average reference
%7) taking off white matter channels
%8) defining trial info
%9) segmenting data
%10) selecting epochs of interest (e.g. stimulus)
%11) exploring signal trialwise to detect artifacts
%12) saving preprocessed clean data

for bb = 1:params.nblocks;
    
    
    %% some parameters
    eval(sprintf('params.fileECOG    = ''block00%d.edf'';', bb))
    eval(sprintf('params.ssIDb        = ''%s_b%d'';', params.ssID, bb)); 
    eval(sprintf('params.tfafile = ''%s_preproTFA_AVG''', params.ssIDb));
        
   
        %% loading data and applying filters
        
        cfg = [];     %initialize configuration file
        cfg.dataset     = [params.pathECOG, params.fileECOG];
        cfg.channel     =['all', params.nonphyschan, params.unknownchan] % non-physiological and unknown channels are removed
        cfg.continuous  = 'yes';
        cfg.dftfreq  = params.notch;
        cfg.bpfreq        = params.bandpass;
        [~,~,ext] = fileparts(params.fileECOG);
        
        data_ecog  = ft_preprocessing(cfg);
        
        %% clean channel names
        
        switch params.cleanchannnames
            case 'yes'
                data_ecog.label = strrep(data_ecog.label, 'EEG ', '');
                data_ecog.hdr.label = strrep(data_ecog.hdr.label,'EEG ', '');
                data_ecog.label = strrep(data_ecog.label, 'MISC ', '');
                data_ecog.label = strrep(data_ecog.label, '_46', '');
                data_ecog.hdr.label = strrep(data_ecog.hdr.label, 'MISC ', '');
            case 'no'
        end
        
        %% loading event channel
        
        cfg = [];
        cfg.dataset     = [params.pathECOG, params.fileECOG];
        cfg.channel     = params.trigchan;        %in this case, cfg is only the trigger channel
        cfg.continuous  = 'yes';
        trig  = ft_preprocessing(cfg);       %raw data of the trigger channel are preprocessed and named trig
        
        %% checking  triggers
        
        triggerdata = trig.trial{1};
        [datPsy,txtPsy] = xlsread([params.pathbehav, params.logfilename]);   %load psychopy data to call the function checking triggers
        [SEEGtriggers, Triggers] = BP_checking_triggers(triggerdata, datPsy, txtPsy, params.fileECOG);
        
        close
        
        %% explore signal for epi artifact
        % channelwise
            for cc = 1: length(data_ecog.label)
                plot(data_ecog.time{1,1}, data_ecog.trial{1,1}(cc,:))
                title(data_ecog.label{cc})
                title(data_ecog.label{cc})
                chan2exclude(cc,1) = str2double(input('Exclude?:', 's')); %1 is for yes
                pause
                
            end
                        
            %% take off bad channels
                
                chan2exclude_lab = data_ecog.label(chan2exclude == 0)';
                
                for ii = 1:length(chan2exclude_lab)
                    str = ['-', chan2exclude_lab{1,ii}];
                    chan2exclude_lab{1,ii} = join(str);
                end
                
                cfg = [];
                cfg.channel     =['all',  chan2exclude_lab] ;
                cfg.continuous  = 'yes';
                data_ecog_clean  = ft_preprocessing(cfg, data_ecog);
            
        %% average reref
        cfg = [];
        cfg.channel = ['all'] ;
        cfg.reref               = 'yes';
        cfg.refchannel          = ['all']; % the average of all the channels uploaded (in the previois cfg.channel) is used ad reference
        
        data_ecogAVG = ft_preprocessing(cfg,data_ecog_clean);         %the data_ecog (see 1st step) are preprocessed with average reference and named data_prepro_erp
        
        %% take WM channels off
        cfg = [];
        cfg.dataset = data_ecogAVG;
        cfg.channel     =['all', params.WMchan];
        data_ecogAVG_clean  = ft_preprocessing(cfg, data_ecogAVG);
        
        %% calling function to define trial info
        
        TrialInfo = importdata([params.pathECOG, params.filePSY]);
        Measures = importdata([params.pathECOG, params.measures]);    %import hand measures
        
        [Trl, TrlCell] = BP_populatetrialinfo(TrialInfo, Measures);           %define trial info (condition, RT, accuracy...)
        [trl] = BP_trialfun_content(cfg,params.fileECOG,trig, SEEGtriggers, Trl, params);                           %define trl (start, duration, end + additional infos)
        

        %% segmenting data
        cfg = [];
        cfg.trl = trl;
        
        data_segm_AVG_all = ft_redefinetrial(cfg, data_ecogAVG_clean);
        timepoint_zero = length(params.epoch(1)*1000:params.epoch(2)*1000) - params.epoch(2)*1000 -1;
        
        %% selecting epochs of interest
        % triggers were set on fixation cross and stimulus onset. I need to
        % select stimulus trials (all the even epochs).
        
        cfg = [];
        cfg.trials = 2:2:length(data_segm_AVG_all.trial);
        data_segm_AVG = ft_selectdata(cfg, data_segm_AVG_all)
        
        
        %% explore trialwise
        
            ntimes = size(data_segm_AVG.time{1,1},2);
            ntrials = length(data_segm_AVG.trial);
            nchan = length(data_segm_AVG.label);
        if cleanresp == 1 | cleanresp == 0.5
            needfigure = str2num(input('Need trial figure?:', 's')); 
            if needfigure == 1
            for cc = 1: nchan
               cc
               f = figure(1)
                set(f, 'Position', [100 -100  2000 2000]);
                % min and max of the non-segmented data to set the axes
                minval = min(data_ecogAVG_clean.trial{1,1}(cc,:));
                maxval =max(data_ecogAVG_clean.trial{1,1}(cc,:));
                for tt = 1:ntrials
                    subplot(7,8,tt);
                    plot(data_segm_AVG.time{1,1}, data_segm_AVG.trial{1,tt}(cc,:))
                    title([ num2str(tt)], 'FontSize', 5 )
                    ylim([minval maxval]);
                    set(gca,'xtick',[]);
                    set(gca,'ytick',[]);
                    
                end
                pause
            end 
            else
            end
            trials2remove = str2num(input('Numbers to exclude?:', 's')); 
            % use selectdata to remove the bad trials 
            cfg = [];
            trialvec = 1:ntrials;
            kkk = 1;
            for ii = 1:ntrials
                if find(trials2remove == trialvec(ii));
                else
                    cfg.trials(1,kkk) = kkk;
                    kkk = kkk+1;
                end
            end
            data_segm_AVG_trclean = ft_selectdata(cfg, data_segm_AVG)
        else
          load([params.OutPath params.ssIDb 'trials2remove']); % if I've already done it, I just upload trials2remove
            cfg = [];
            trialvec = 1:ntrials;
            kkk = 1;
            for ii = 1:ntrials
                if find(trials2remove == trialvec(ii));
                else
                    cfg.trials(1,kkk) = kkk;
                    kkk = kkk+1;
                end
            end
            data_segm_AVG_trclean = ft_selectdata(cfg, data_segm_AVG)
        end
    
    save([params.OutPath params.ssIDb 'preproc_TFA_trclean'], 'data_segm_AVG_trclean');
    if cleanresp == 1 | cleanresp == 0.5
        exist trials2remove
        if ans == 1
            save([params.OutPath params.ssIDb 'trials2remove'],'trials2remove');
        end 
    else
    end
end


save([params.OutPath params.ssIDb 'preproc_TFA_clean'], 'data_segm_AVG_trchclean');

















