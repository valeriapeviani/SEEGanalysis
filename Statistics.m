%% STATISTICS
% apply permutation statistics on high-gamma and beta (or whatever you want) power.
% 1st step --> select responsive channels (channels in which "task"
% activity is significantly higher or lower than "baseline activity"
 


%% some more parameters
    
    % convert baseline time into indicesf
    timepoint_zero = length(params.epoch(1)*1000:params.epoch(2)*1000) - params.epoch(2)*1000 -1;
    baseidx = [timepoint_zero + baseline_window(1) timepoint_zero + baseline_window(2)];
    taskidx =  [timepoint_zero + task_window(1) timepoint_zero + task_window(2)];
    
    for cond = 1: length(params.conditions)
    % load matrix: nchan x frequency range  x  trials x times
        
        load([params.OutPath 'hilbTFdB_' params.conditions{cond}]);
        eval(sprintf('ntrials%s = size(tfhilb_DB_%s, 3);',params.conditions{cond}, params.conditions{cond}));
        eval(sprintf('nchans = size(tfhilb_DB_%s, 1);', params.conditions{cond}));
        eval(sprintf('ntimes = size(tfhilb_DB_%s, 4);', params.conditions{cond}));
    end  
        
     
    
    % load 1 block just for the channel labels 

for bb = 1; 
eval(sprintf('block%d = load([params.OutPath params.ssID ''_b'' num2str(bb) ''preproc_TFA_clean'']);',  bb));
end  


%% compare baseline (pres-timulus) vs task (post-stimulus) to select "responsive channels for further analysis

%% prepare data
for ff = 1:length(params.targetfrex)
    test_table = zeros(nchans, 8);
    for cc = 1:nchans
        for cond = 1: length(params.conditions)
            % create temporary matrix for the two conditions
            if cond ==1 ;
                datatemp1 = zeros(ntrialshand, ntimes);
                    eval(sprintf('datatemp1 = squeeze(tfhilb_DB_%s(cc, ff, :, :));', params.conditions{cond}));
            else
                datatemp2 = zeros(ntrialsmob, ntimes);
                    eval(sprintf('datatemp2 = squeeze(tfhilb_DB_%s(cc, ff, :, :));', params.conditions{cond}));
            end
        end
        
        data_all = zeros(ntrialshand + ntrialsmob, ntimes);
        data_all = cat(1,datatemp1, datatemp2);
        
        % prepare matrices to compare task vs baseline
        data_all_bas = zeros(ntrialshand + ntrialsmob,1);
        data_all_task = zeros(ntrialshand + ntrialsmob,1);
        data_all_bas = squeeze(mean(data_all(:, baseidx(1):baseidx(2)),2));
        data_all_task = squeeze(mean(data_all(:, taskidx(1):taskidx(2)),2));
        
        
        %% permutation statistics
        for cc = 1:nchans
            for cond = 1: length(params.conditions)
                % create temporary matrix for the two conditions
                if cond ==1 ;
                    datatemp1 = zeros(ntrialshand, ntimes);
                        eval(sprintf('datatemp1 = squeeze(tfhilb_DB_%s(cc, ff, :, :));', params.conditions{cond}));
                else
                    datatemp2 = zeros(ntrialsmob, ntimes);
                        eval(sprintf('datatemp2 = squeeze(tfhilb_DB_%s(cc, ff, :, :));', params.conditions{cond}));
                end
            end
            
            data_all = zeros(ntrialshand + ntrialsmob, ntimes);
            data_all = cat(1,datatemp1, datatemp2);
            
            % prepare matrices to compare task vs baseline
            data_all_bas = zeros(ntrialshand + ntrialsmob,1);
            data_all_task = zeros(ntrialshand + ntrialsmob,1);
            data_all_bas = squeeze(mean(data_all(:, baseidx(1):baseidx(2)),2));
            data_all_task = squeeze(mean(data_all(:, taskidx(1):taskidx(2)),2));
            % compute true difference between baseline and task
            true_bastask_diff = mean(data_all_task) - mean(data_all_bas);
            
            % pool data together
            data_bastask = zeros(1, (ntrialshand + ntrialsmob) * 2);
            data_bastask = cat(1, data_all_bas, data_all_task);
            % truelabels = 1 = baseline, 2 = task
            truelabels_bastask = zeros(1, (ntrialshand + ntrialsmob) * 2);
            truelabels_bastask = [repmat(1, ntrialshand + ntrialsmob, 1); repmat(2, ntrialshand + ntrialsmob, 1)];
            niterations = 1000;
            % % compute true difference between baseline and task (TASK - BASELINE)
            %
            for permi = 1:niterations
                
                % shuffle the labels (0 = baseline, 1 = task)
                shuflabels_bastask = truelabels_bastask(randperm(length(truelabels_bastask)));
                
                % compute the mean difference of the shuffled labels
                permshuf_bastask_diff(permi,1) = mean(data_bastask(shuflabels_bastask==2)) - mean(data_bastask(shuflabels_bastask==1));
                
            end
            
            %                         %show the distribution (to check for normality)
            %                         figure(1), clf, hold on
            %                         histogram(permshuf_bastask_diff,40)
            %                         %plot([1 1]*true_bastask_diff,get(gca,'ylim'),'r--','linew',3)
            %                         plot([1 1]*0.3,get(gca,'ylim'),'r--','linew',3)
            %                         legend({'Shuffled labels';'Observed'},'Location', 'northwest')
            %                         set(gca,'xlim',[-1 1])
            %                         xlabel('Difference under H_0')
            %                         ylabel('Count')
            
            %
            % compute pvalue (how much the true value diverges from the H0 distribution
            % (obtained by shuffling the labels)
            permmean = mean(permshuf_bastask_diff);
            permstd  = std(permshuf_bastask_diff);
            
            % formula for z-score
            zdist(cc,:) = (true_bastask_diff - permmean) / permstd ;
            
            % can convert to p-value
            pvals(cc,1) = normpdf(zdist(cc,:));
            
            if pvals(cc,1) <.05
                hperm(cc,1) = 1;
            else
                hperm(cc,1) = 0;
            end
        end
        test_tableT = array2table(test_table,'VariableNames',{'perm_p', 'permZ', 'permHP', 'adjpermp', 'adjpermHP', 'chanlabels', 'fdrthresh_ttest', 'fdrthresh_perm'});
        test_tableT.perm_p(:) = pvals;
        test_tableT.permZ(:) = zdist;
        test_tableT.permHP(:) = hperm;
        
        for corri = 1:length(params.cor1threshperm)
            
            [h, ~, ~, adj_p]=fdr_bh(pvals, params.cor1threshperm(corri));
            test_tableT.adjpermp(:)= adj_p;
            test_tableT.adjpermHP(:)=  h;
            test_tableT.fdrthresh_perm(:) = params.cor1threshperm;
        end
        
        %                     % save output
            writetable(test_tableT, [params.OutPath params.ssID 'Stats_HILB_BasVSTask_fdr' num2str(params.cor1thresh(corri)*1000) '_' num2str(params.cor1threshperm(corri)*1000) '_' params.targetfrex{ff} '.xls']);
        %% selecting permutation-based activated channels
            [ncontacts, txtcontacts, respcontacts] = xlsread([params.OutPath params.ssID 'Stats_HILB_BasVSTask_fdr' num2str(params.cor1thresh(corri)*1000) '_' num2str(params.cor1threshperm(corri)*1000) '_' params.targetfrex{ff} '.xls'],1);
        
        load([params.pathECOG, params.blockord]);   %load block order
        
        
        %% select channels based on permutation statistics
        statcrit = 2;
        HPidx = 4;
        
        ncont = nchans;
        cc = 1;
        for ch = 1:ncont;
            if ncontacts(ch,HPidx) ==1;
                eval(sprintf('respcont_%s(cc,1) = chanlabels(ch,1);', params.targetfrex{ff}));
                eval(sprintf('respcontIND_%s(cc,1) = find(strcmp(chanlabels, respcont_%s(cc,1)));', params.targetfrex{ff},params.targetfrex{ff}));
                cc = cc+1;
                
            end
        end
        
        eval(sprintf('exist respcont_%s;', params.targetfrex{ff}));
        if ans == 1
            
            
            eval(sprintf('nactchans = size(respcont_%s,1);', params.targetfrex{ff}));
            
            %% save resp contacts
            if statcrit == 1
                eval(sprintf('writetable(cell2table(respcont_%s), [params.OutPath params.ssID ''_sel_cont_ttest_handmob'' params.targetfrex{ff}], ''FileType'', ''spreadsheet'');',  params.targetfrex{ff}));
            else
                eval(sprintf('writetable(cell2table(respcont_%s), [params.OutPath params.ssID ''_sel_cont_permut_handmob'' params.targetfrex{ff}], ''FileType'', ''spreadsheet'');',  params.targetfrex{ff}));
            end
            
            handsel = zeros(ntrialshand, 1);
            mobsel = zeros(ntrialsmob, 1);
            zdist_handmob = zeros(nactchans,1);
            pvals_handmob = zeros(nactchans,1);
            h_handmob = zeros(nactchans,1);
            adjp_handmob = zeros(nactchans,1);
            
            % prepare matrices and average over time
            for acc = 1: nactchans
                if wav_or_hilb == 1
                    eval(sprintf('handsel(:,1) = squeeze(mean(tf_DB_trials_hand(respcontIND_%s(acc), ff, :, timepoint_zero:timepoint_zero+task_window(2)),4));', params.targetfrex{ff}));
                    eval(sprintf('mobsel(:,1) = squeeze(mean(tf_DB_trials_mob(respcontIND_%s(acc), ff, :, timepoint_zero:timepoint_zero+task_window(2)),4));', params.targetfrex{ff}));
                else
                    eval(sprintf('handsel(:,1) = squeeze(mean(tfhilb_DB_hand(respcontIND_%s(acc), ff, :, timepoint_zero:timepoint_zero+task_window(2)),4));', params.targetfrex{ff}));
                    eval(sprintf('mobsel(:,1) = squeeze(mean(tfhilb_DB_mob(respcontIND_%s(acc), ff, :, timepoint_zero:timepoint_zero+task_window(2)),4));', params.targetfrex{ff}));
                end
                
                %% run permutation statistics
                
                condmat_handmob = [zeros(ntrialshand,1); repmat(1,ntrialsmob,1)];
                % 0 = hand, 1 = mob
                niterations = 1000;
                
                % compute true difference between conditions
                true_handmob_diff = mean(handsel) - mean(mobsel);
                % pool the data
                datapooled = [handsel; mobsel];
                
                for permi = 1:niterations
                    
                    % shuffle the labels (0 = hand, 1 = mobile)
                    shufcondmat_handmob = condmat_handmob(randperm(length(condmat_handmob)));
                    
                    % compute the mean difference of the shuffled labels
                    permshuf_handmob_diff(permi,1) = mean(datapooled(shufcondmat_handmob==0)) - mean(datapooled(shufcondmat_handmob==1));
                    
                end
                
                % compute pvalue (how much the true value diverges from the H0 distribution
                % (obtained by shuffling the labels)
                permmean = mean(permshuf_handmob_diff);
                permstd  = std(permshuf_handmob_diff);
                
                % formula for z-score
                zdist_handmob(acc,1) = (true_handmob_diff - permmean) / permstd ;
                
                % can convert to p-value
                pvals_handmob(acc,1) = normpdf(zdist_handmob(acc,1));
                
            end
            
            % put results in a table
            test_handmob = zeros(nactchans, 5);
            test_handmob(:,1) = zdist_handmob;
            test_handmob(:,2) = pvals_handmob;
            
            [h_handmob, ~, ~, adjp_handmob]=fdr_bh(pvals_handmob(:), 0.05);
            
            test_handmob(:,4) = h_handmob;
            test_handmob(:,3) = adjp_handmob;
            
            test_handmobT = array2table(test_handmob, 'VariableNames', {'zscore', 'pvalue', 'adjp', 'adjH', 'chan'});
            eval(sprintf('test_handmobT.chan = respcont_%s;', params.targetfrex{ff}));
            
            %% save data
                writetable(test_handmobT, [params.OutPath params.ssID 'conditions_compar_HILB_permut_' params.targetfrex{ff} '.xls']);
            
        else
        end
    end
    
    end 

