%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract high-gamma and beta through convolution (complex morelet wavelet) 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) load .mat files
% 2) define TFA paramteres
% 3) run TFA via complex morelet wavelets in the frequency domain 
% to extract 1) a matrix nchannels x nconditions  x  nfrequencies x ntimes (for each
% freq band (e.g. high-gamma)) in which the trials are averaged to plot the
% freq spectra, 2) a matrix nchannels x frequencyrange (e.g. high-gamma)  x ntrials x ntimes 
% for each condition, in which the power values for each frequency band are
% averaged together, to plot average high-value over time. 
% 4) power normalization 
% 5) frequency spectra for each condition and frequency range
% 6) plots of the power of the two conditions averaged over trials and
% frequencies 
% 7) save data

%% load preprocessed files

for bb = 1:params.nblocks    
eval(sprintf('block%d = load([params.OutPath params.ssID ''_b'' num2str(params.nblocks(bb)) ''preproc_TFA_allclean'']);',  params.nblocks(bb)));
eval(sprintf('block%d = block%d.data_segm_AVG_trchclean;', params.nblocks(bb),params.nblocks(bb)));
end

%% general TFA parameters
% looping over the two frequency bands I'm interested in 
for ff = 1:length(params.targetfrex)

% prepare for convolution in the frequency domain
eval(sprintf('min_frex = params.frex_%s(1);', params.targetfrex{ff}));
eval(sprintf('max_frex = params.frex_%s(3);', params.targetfrex{ff}));
eval(sprintf('num_frex = params.frex_%s(2);', params.targetfrex{ff}));
range_cycles = params.range_cycles;
frex = linspace(min_frex, max_frex, num_frex);% range of wavelet cycles
nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex); % n of cycles in logspace
baseline_window = [ -500 -200 ];
% convert baseline time into indices
timepoint_zero = length(params.epoch(1)*1000:params.epoch(2)*1000) - params.epoch(2)*1000 -1;
%tpz = length(-1000:1500) - 1500 
baseidx = [timepoint_zero + baseline_window(1) timepoint_zero + baseline_window(2)];

time = params.epoch(1):1/params.srate:params.epoch(2);
half_wave = (length(time)-1)/2;
ntimes = size(block1.time{1,1},2);

%% pool together data for the two conditions
nchans = length(block1.label);  % number of channels does not change across blocks 
    % loop over channels
for cc = 1:nchans
    
    % put data in matrices
    for bb = params.nblocks'
        eval(sprintf('ntrials = length(block%d.trial);', params.nblocks(bb)));
        for nt = 1:ntrials
            eval(sprintf('block%d_dat(nt,:) = block%d.trial{1,nt}(cc,:) ;', params.nblocks(bb), params.nblocks(bb)));
        end
    end

    % loop over conditions
for cond = 1:length(params.conditions);
    
    % compute number of trials per condition
    if cond == 1
        eval(sprintf('ntrials(cond) = length(block%d.trial) + length(block%d.trial);', params.nblocks(1), params.nblocks(2)));
    elseif cond == 2
        eval(sprintf('ntrials(cond) = length(block%d.trial) + length(block%d.trial);', params.nblocks(3), params.nblocks(4)));
    end
        alldata = zeros(ntrials(cond), ntimes);

    % create a matrix for each condition (each block, each channel)
    if cond == 1
        eval(sprintf('alldata(:,:) = [block%d_dat; block%d_dat] ;', params.nblocks(1), params.nblocks(2)));
    elseif cond == 2
        eval(sprintf('alldata(:,:) = [block%d_dat; block%d_dat] ;', params.nblocks(3), params.nblocks(4)));
    end
    
        
%% channel-specific FFT parameters
nKern = length(time);
nData = length(time)*ntrials(cond); %because I'm gonna reshape it in a continuous signal
nConv = nKern+nData-1;

% initialize output time-frequency data
tf = zeros(length(frex),length(time));

%reshape data matrix in concatenating all the trials in a row 
data_resh = reshape(alldata(:,:),1,[]);

% FFT of signal (doesn't change on frequency iteration)
dataX = fft(data_resh ,nConv);  % nConv specifies the zero padding

% loop over frequencies  
as_resh_allfr = zeros(length(frex), ntrials(cond), ntimes);

    for fi=1:length(frex)
        % compute parameter s for the gaussian
        s = nCycles(fi)/(2*pi*frex(fi));
        % compute complex morelet wavelet (complex sine wave .* gaussian)
        cmw  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
        % fft of the kernel (complex morelet wavelet)
        cmwX = fft(cmw,nConv);
        % demaximize spectrum 
        cmwX = cmwX./max(cmwX);
        
        % run convolution, trim edges, and reshape to 2D (time X trials)
        as = ifft(cmwX.*dataX);
        as_cut = as(half_wave:end-(half_wave+1));
        as_resh = reshape(as_cut,ntrials(cond),ntimes);
        
        % average power data ACROSS TRIALS and put it into big matrix
        tf_freq(fi,:) = mean(abs(as_resh).^2,1);
        % cumulate frequencies
        as_resh_allfr(fi,:,:) = as_resh;
    end
    
    %% power normalization
    % for the matrix averaged over trials (keeping frequencies) or for the
    % matrix that will be averaged over frequencies (keeping trials)
    matrices = {'freq', 'trials'};
    
    for ii = 1:length(matrices)   % two different matrices: tf is averaged over trials, tf_trials is averaged over frex 
    if ii == 1;
        % power normalization when the data are already averaged over
        % trials
    activity = tf_freq;
    baseline = mean(tf_freq(:,baseidx(1):baseidx(2)),2);
    tfDB = 10*log10( bsxfun(@rdivide, activity, baseline));
    tfpc = bsxfun(@rdivide, (100.*(activity - baseline)), baseline) ;
    cutidx = 50; % to remove edge artifacts
    tf2plot = zeros(size(activity,1), size(activity(1,1:end-cutidx),2));
    tfDB2plot = zeros(size(activity,1), size(activity(1,1:end-cutidx),2));
    tfpc2plot = zeros(size(activity,1), size(activity(1,1:end-cutidx),2));
    tf2plot(:,:) = activity(:,1:end-cutidx);
    tfDB2plot(:,:) = tfDB(:,1:end-cutidx);
    tfpc2plot (:,:) = tfpc(:,1:end-cutidx);
    
    else
        % power normalization within reach trial (afterwards, I can average
        % the frequency together)
    for iii = 1:ntrials(cond);
    activity(:,iii,:) = abs(as_resh_allfr(:,iii,:)).^2;
    % baseline is computed trialwise
    baselineeachtr = zeros(length(frex),ntrials(cond));
    baselineeachtr(:,iii) = mean(activity(:,iii,baseidx(1):baseidx(2)),3);
    tfDB(:,iii,:) = 10*log10( bsxfun(@rdivide, activity(:,iii,:), baselineeachtr(:,iii)));
    tfpc(:,iii,:) = bsxfun(@rdivide, (100.*(activity(:,iii,:) - baselineeachtr(:,iii))), baselineeachtr(:,iii)) ;
    cutidx = params.cutidx; % to remove edge artifacts
    tf2plot = zeros(size(activity,1), ntrials(cond), size(activity(1,iii,1:end-cutidx),3));
    tfDB2plot = zeros(size(activity,1),ntrials(cond), size(activity(1,iii,1:end-cutidx),3));
    tfpc2plot = zeros(size(activity,1),ntrials(cond), size(activity(1,iii,1:end-cutidx),3));
    % cumulate results recovering trials
    tf2plot(:,iii,:) = activity(:,iii,1:end-cutidx);
    tfDB2plot(:,iii,:) = tfDB(:,iii,1:end-cutidx);
    tfpc2plot (:,iii,:) = tfpc(:,iii,1:end-cutidx);    
    end
    
    % average power data ACROSS FREQUENCIES
    tf_trials = zeros(ntrials(cond), ntimes);
    tf_trials(:,:) = mean(abs(as_resh_allfr).^2,1);
    tfDB_trials = zeros(ntrials(cond), ntimes);
    tfDB_trials(:,:) = mean(tfDB,1);
    tfPC_trials = zeros(ntrials(cond), ntimes);
    tfPC_trials(:,:) = mean(tfpc,1);
    end 
    
    %% create files to be saved later
    
    if ii == 1
        % matrix nchan x condition  x  n of frequencies x times (for each
        % freq range) AVERAGED OVER TRIALS
    eval(sprintf('tf_raw_%s_%s(cc,cond,:,:) = tf_%s(:,:); ', matrices{ii}, params.targetfrex{ff}, matrices{ii}));
    eval(sprintf('tf_DB_%s_%s(cc,cond,:,:) = tfDB(:,:); ', matrices{ii}, params.targetfrex{ff}));
    eval(sprintf('tf_PC_%s_%s(cc,cond,:,:) = tfpc(:,:); ', matrices{ii}, params.targetfrex{ff}));
    else
        % matrix nchan x frequency range  x  trials x times (for each freq
        % range) AVERAGED OVER FREQUENCIES
    eval(sprintf('tf_raw_%s_%s(cc,ff,:,:) = tf_%s(:,:); ', matrices{ii}, params.conditions{cond}, matrices{ii}));
    eval(sprintf('tf_DB_%s_%s(cc,ff,:,:) = tfDB_trials(:,:); ', matrices{ii}, params.conditions{cond}));
    eval(sprintf('tf_PC_%s_%s(cc,ff,:,:) = tfPC_trials(:,:); ', matrices{ii}, params.conditions{cond}));
    end 
    
    
    %% spectra, averaged over trials with raw signal, DB corrected and percent change corrected 

    if ii == 1
scrsz = get(groot,'ScreenSize');

figure(1), clf
set(gcf, 'Position',[1 1 scrsz(3)/2 scrsz(4)])

subplot(3,1,1)
contourf(time(1:end-cutidx),frex,tf2plot,40,'linecolor','none')
hold on
title(['Raw ' params.targetfrex{ff} 'spectra chan ' block1.label{cc} 'cond ' params.conditions{cond}]);
xlabel('Time (ms)'), 
xticks(time(1):0.2:time(end-cutidx));
xticklabels(time(1)*1000:200:time(end-cutidx)*1000);
plot([0 0 ], [min_frex max_frex], '--r'); 
ylabel('Frequency (Hz)')
if ff == 1
yticks(min_frex:2:max_frex);
else
yticks(min_frex:(max_frex-min_frex)/10:max_frex);
end 
c = colorbar;
c.Label.String = 'Raw Power (%change)'; 


subplot(3,1,2)
contourf(time(1:end-cutidx),frex,tfDB2plot,40,'linecolor','none')
hold on
title(['DBnorm ' params.targetfrex{ff} ' spectra chan ' block1.label{cc} 'cond ' params.conditions{cond}]);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
xticks(time(1):0.2:time(end-cutidx));
xticklabels(time(1)*1000:200:time(end-cutidx)*1000);
plot([0 0 ], [min_frex max_frex], '--r'); 
if ff == 1
yticks(min_frex:2:max_frex);
else
yticks(min_frex:(max_frex-min_frex)/10:max_frex);
end 
c =colorbar;
c.Label.String = 'Normalized Power (%change)'; 

subplot(3,1,3)
contourf(time(1:end-cutidx),frex,tfpc2plot,40,'linecolor','none')
hold on
title(['PCnorm ' params.targetfrex{ff} ' spectra chan ' block1.label{cc} 'cond ' params.conditions{cond}]);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
xticks(time(1):0.2:time(end-cutidx));
xticklabels(time(1)*1000:200:time(end-cutidx)*1000);
plot([0 0 ], [min_frex max_frex], '--r'); 
if ff == 1
yticks(min_frex:2:max_frex);
else
yticks(min_frex:(max_frex-min_frex)/10:max_frex);
end 
c = colorbar;
c.Label.String= 'Normalized Power  (dB)'; 

saveas(gcf, [params.outpath_plots '\' params.targetfrex{ff} '\threespectra' '\TFA' params.targetfrex{ff} '_' block1.label{cc} '_' params.conditions{cond}], 'jpg');


%% spectra with only normalized power (dB)

figure(1), clf
set(gcf, 'Position',[1 1 scrsz(3) scrsz(4)/3])
contourf(time(1:end-cutidx),frex,tfDB2plot,40,'linecolor','none')
hold on
title(['DBnorm ' params.targetfrex{ff} ' spectra chan ' block1.label{cc} 'cond ' params.conditions{cond}]);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
xticks(time(1):0.2:time(end-cutidx));
xticklabels(time(1)*1000:200:time(end-cutidx)*1000);
plot([0 0 ], [min_frex max_frex], '--r'); 
if ff == 1
yticks(min_frex:2:max_frex);
else
yticks(min_frex:(max_frex-min_frex)/10:max_frex);
end 
c =colorbar;
c.Label.String = 'Normalized Power (%change)'; 

saveas(gcf, [params.outpath_plots '\' params.targetfrex{ff} '\onlyDB' '\TFA' params.targetfrex{ff} '_' block1.label{cc} '_' params.conditions{cond}], 'jpg');
close 

    else
    end  
clear  activity baseline tfDB tfpc tf2plot tfDB2plot tfpc2plot
    
    end 
end 

%% plot power averaged over trials and frequencies 
figure(2), clf
options.handle     = figure(2);

for cond = 1:length(params.conditions);
eval(sprintf('data2plot_%s = zeros(ntrials(cond), ntimes);', params.conditions{cond}));
eval(sprintf('data2plot_%s(:, :)  = squeeze(tf_DB_trials_%s(cc,ff,:,:));',params.conditions{cond}, params.conditions{cond}));

% average over trials
eval(sprintf('data2plot_av_%s = mean(data2plot_%s,1);', params.conditions{cond}, params.conditions{cond}));
eval(sprintf('data2plot_std_%s = std(data2plot_%s,1)'';', params.conditions{cond}, params.conditions{cond}));
% smooth 
sparam = params.smoothpar(ff);
eval(sprintf('data2plot_av_smooth_%s = smooth(data2plot_av_%s, sparam, ''moving'');', params.conditions{cond}, params.conditions{cond}));
eval(sprintf('data2plot_std_smooth_%s = smooth(data2plot_std_%s, sparam, ''moving'');', params.conditions{cond}, params.conditions{cond}));

plot_tw = [timepoint_zero - 600 timepoint_zero + 800];  % specify temporal window to be plotted
%xval = [0:100:1200];
% if ff == 1 % beta
%     ylimits = [0.9 1.5];
%     xlimits = [0 1200];
% elseif ff ==2 % HG
%     ylimits = [0.9 1.3];
%     xlimits = [0 1200];
% end 
colors_area = [[102 153 255]./255; [255 53 0]./255]; % red and blue
colors_line = [[51 0 255]./255; [204 0 0 ]./255];
% options for shaded area error plot
options.alpha      = 0.2;
options.line_width = 1;
options.error      = 'std'; % i want to plot standard deviation
options.ntrials     = ntrials(cond);
options.color_area = colors_area(cond,:);
options.color_line = colors_line(cond,:);

eval(sprintf('plot_areaerrorbar(data2plot_av_smooth_%s(plot_tw(1):plot_tw(2))'', data2plot_std_smooth_%s(plot_tw(1):plot_tw(2))'', options );', params.conditions{cond}, params.conditions{cond}));
hold on
% plot timewindow lines
xst = [400 400]; % vertical line to signal time 0
xnd = [1000 1000 ]; % vertical line to signal end of time window considered in the analysis (600 ms after the stimulus onset)
y = [-2 2];
plot(xst, y,'LineWidth', 1,'LineStyle', '--', 'Color', 'r', 'HandleVisibility','off');
plot(xnd, y,'LineWidth', 1,'LineStyle', '--', 'Color', 'r', 'HandleVisibility','off');
end 
% other options
xlim([0 plot_tw(2)-plot_tw(1)]);
xticklabels([600 - timepoint_zero:200:timepoint_zero]);
xlabel('msec');
%ylim([-2 2]);
[leg,labelhandles,outH,outM] =legend('Hand','Mobile');%, 'Location')%, legendpos);
legendpos = 'northeast';
title(['average ' params.targetfrex{ff} ' power ' block1.label{cc}]);
ylabel(['normalized ' params.targetfrex{ff} ' power (dB)']);
    
saveas(gca, [params.outpath_plots '/' params.targetfrex{ff} '/averagepower/'  params.ssID '_Average_' params.targetfrex{ff} '_power_' block1.label{cc}], 'png' );
%close 

 
end 


end 
%% save data
for ii = 1:length(matrices)
    if ii == 1;
        for ff = 1:length(params.targetfrex);
        % matrix nchan x condition  x  n of frequencies x times (for each
        % freq range) AVERAGED OVER TRIALS
    eval(sprintf('save([params.outpath_plots ''TFraw_avovertrials_'' params.targetfrex{ff}], ''tf_raw_%s_%s''); ', matrices{ii}, params.targetfrex{ff}));
    eval(sprintf('save([params.outpath_plots ''TFdB_avovertrials_'' params.targetfrex{ff}], ''tf_DB_%s_%s''); ', matrices{ii}, params.targetfrex{ff}));
    eval(sprintf('save([params.outpath_plots ''TFpc_avovertrials_'' params.targetfrex{ff}], ''tf_PC_%s_%s''); ', matrices{ii}, params.targetfrex{ff}));
        end
    else 
        % matrix nchan x frequency range  x  trials x times (for each freq
        % range) AVERAGED OVER FREQUENCIES
        for cond = 1:length(params.conditions);
    eval(sprintf('save([params.outpath_plots ''TFraw_avoverfrex_'' params.conditions{cond}], ''tf_raw_%s_%s''); ', matrices{ii}, params.conditions{cond}));
    eval(sprintf('save([params.outpath_plots ''TFdB_avoverfrex_'' params.conditions{cond}], ''tf_DB_%s_%s''); ', matrices{ii}, params.conditions{cond}));
    eval(sprintf('save([params.outpath_plots ''TFpc_avoverfrex_'' params.conditions{cond}], ''tf_PC_%s_%s''); ', matrices{ii}, params.conditions{cond}));
        end
    end 
end 



