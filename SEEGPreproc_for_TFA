%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing data for TFA
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) loading EDF data without non-phys channels and unknown channels
% 1a) load trigger channel, check sync between triggers and logfile
% timestamps --> here custom-made function "Check_triggers" is needed
% 2) cleaning channelwise (optional) 
% 3) average reference
% 4) remove white matter channels
% 5) segment data --> here custom-made functions "Populatetrialinfo" and
% "Trialfun_content" are needed
% 6) cleaning trialwise 
% 7) cleaning channelwise (again) 
%cleaning trialwise, cleaning channelwise again (to remove channels with many noisy trials)
% save preprocessed data in .mat format 

%%
%% upload blocks 
cleanresp = str2num(input('Need to clean?:', 's'));
% 0 = no need to clean anything (or already cleaned). It goes to channel
% rejection after the trials. 
% 0.5 = clean only trialwise
% 1 = clean also the channelwise (before cleaning the trials)

for bb = 1:params.nblocks;
    
    
    %% some parameters
    eval(sprintf('params.fileECOG    = ''block00%d.edf'';', bb))
    eval(sprintf('params.ssIDb        = ''%s_b%d'';', params.ssID, bb)); 
    eval(sprintf('params.tfafile = ''%s_preproTFA_AVG''', params.ssIDb));
       
        %% loading data and applying filters
        
        cfg = [];     %initialize configuration file
        cfg.dataset     = [params.pathSEEG, params.fileECOG];
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
        cfg.dataset     = [params.pathSEEG, params.fileECOG];
        cfg.channel     = params.trigchan;        %in this case, cfg is only the trigger channel
        cfg.continuous  = 'yes';
        trig  = ft_preprocessing(cfg);       %raw data of the trigger channel are preprocessed and named trig
        
        %% checking  triggers
        
        triggerdata = trig.trial{1};
        [datPsy,txtPsy] = xlsread([params.pathbehav, params.logfilename]);   %load logfile to call the function checking triggers
        [SEEGtriggers, Triggers] = Check_triggers(triggerdata, datPsy, txtPsy, params.fileECOG);
        
        close
        
        %% explore signal for epi artifact
        % channelwise
        if cleanresp == 1
            for cc = 1: length(data_ecog.label)
                plot(data_ecog.time{1,1}, data_ecog.trial{1,1}(cc,:))
                title(data_ecog.label{cc})
                chan2exclude(cc,1) = str2double(input('Exclude?:', 's')); %1 is for yes
                pause
                
            end
            
            badch = str2double(input('Need to clean?:', 's')); %0 is for yes
            if badch ==1
                disp('Repeat the process above uncommenting chan2exclude');
            end
            
            %% take off bad channels
            if badch == 1
                
                %chan2exclude = repmat(0,1,184);
                %chan2exclude = randi([0 1],1, 184);
                chan2exclude_lab = data_ecog.label(chan2exclude == 0)';
                
                for ii = 1:length(chan2exclude_lab)
                    str = ['-', chan2exclude_lab{1,ii}];
                    chan2exclude_lab{1,ii} = join(str);
                end
                
                cfg = [];
                cfg.channel     =['all',  chan2exclude_lab] ;
                cfg.continuous  = 'yes';
                data_ecog_clean  = ft_preprocessing(cfg, data_ecog);
            else
                data_ecog_clean = data_ecog
            end
            
        elseif cleanresp == 0.5 |cleanresp == 0 ;
            data_ecog_clean = data_ecog   % if it is specified that I don't need 2 clean
        end
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
        
        TrialInfo = importdata([params.pathSEEG, params.filePSY]);
        Measures = importdata([params.pathSEEG, params.measures]);    %import hand measures
        
        [Trl, TrlCell] = Populatetrialinfo(TrialInfo, Measures);           %define trial info (condition, RT, accuracy...)
        [trl] = Trialfun_content(cfg,params.fileECOG,trig, SEEGtriggers, Trl, params);                           %define trl (start, duration, end + additional infos)
        
        %% segmenting data
        cfg = [];
        cfg.trl = trl;
        
        data_segm_AVG_all = ft_redefinetrial(cfg, data_ecogAVG_clean);
        timepoint_zero = length(params.epoch(1)*1000:params.epoch(2)*1000) - params.epoch(2)*1000 -1;
        
        %% selecting stimulus trials
        
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
                    %             hold on
                    %             plot(data_segm_AVG.time{1,1},repmat(meantrial,1, ntimes));
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
    else
    end 
end


%% channel removal (after having looked at all the trials)
chan2exclude_end = str2num(input('Chan indeces to exclude?:', 's')); %1 is for yes
chan2exclude_end_labs = data_segm_AVG_trclean.label(chan2exclude_end);

for bb = block2keep';
    eval(sprintf('params.ssIDb = ''%s_b%d'';', params.ssID, bb)); 
    eval(sprintf('params.tfafile = ''%spreproc_TFA_trclean''', params.ssIDb));
    load([params.OutPath, params.tfafile]);   %load block order
    
    for ii = 1:length(chan2exclude_end_labs)
        chan2exclude_end_labsM{ii,1} = ['-', chan2exclude_end_labs{ii,1}];
        %str = ['-', chan2exclude_end_labs{ii,1}];
        %chan2exclude_end_labs{ii,1} = join(str);
    end
    
        cfg = [];
        cfg.channel     =['all',  chan2exclude_end_labsM'] ;
        cfg.continuous  = 'yes';
        data_segm_AVG_trchclean  = ft_preprocessing(cfg, data_segm_AVG_trclean);
        save([params.OutPath params.ssIDb 'preproc_TFA_allclean'], 'data_segm_AVG_trchclean');

end 

exist chan2exclude_end;
if ans == 1
    save([params.OutPath params.ssID 'chan2exclude_aftertrialrej'],'chan2exclude_end');
end 

save([params.OutPath params.ssID 'block2keep'], 'block2keep');
















