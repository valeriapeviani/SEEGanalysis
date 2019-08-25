%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract high-gamma and beta through filter-hilbert 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load preprocessed data
for bb = params.nblocks 
eval(sprintf('block%d = load([params.OutPath params.ssID ''_b'' num2str(bb) ''preproc_TFA_clean'']);',  bb));
end

%% general TFA parameters
% looping over the two frequency bands I'm interested in 
for ff = 1:length(params.targetfrex)

% built fir filter for hilbert transform 
% filter parameters
srate   = params.srate; % hz
nyquist = srate/2;
frange  = params.frexhilb (ff,:);
transw  = .1;
order   = round( 16*srate/frange(1) );
shape   = [ 0 0 1 1 0 0 ];
frex    = [ 0 frange(1)-frange(1)*transw frange frange(2)+frange(2)*transw nyquist ] / nyquist;

% filter kernel
filtkern = firls(order,frex,shape);
% compute the power spectrum of the filter kernel
filtpow = abs(fft(filtkern)).^2;
% compute the frequencies vector and remove negative frequencies
hz      = linspace(0,srate/2,floor(length(filtkern)/2) + 1);
filtpow = filtpow(1:length(hz));

% % plot the filter kernel
% figure(2), clf
% subplot(131)
% plot(filtkern,'linew',2)
% xlabel('Time points')
% title('Filter kernel (firls)')
% axis square
% 
% % % make the plot look nicer
% %set(gca,'xlim',[0 20])
% xlabel('Frequency (Hz)'), ylabel('Filter gain')
% legend({'Actual';'Ideal'})
% title('Frequency response of filter (firls)')
% axis square
% % 
% subplot(133), hold on
% plot(hz,10*log10(filtpow),'ks-','linew',2,'markersize',10,'markerfacecolor','w')
% plot([1 1]*frange(1),get(gca,'ylim'),'k:')
% set(gca,'xlim',[0 frange(1)*4],'ylim',[-50 2])
% xlabel('Frequency (Hz)'), ylabel('Filter gain (dB)')
% title('Frequency response of filter (firls)')
% axis square

baseline_window = [ -500 -100 ];
% convert baseline time into indices
timepoint_zero = length(params.epoch(1)*1000:params.epoch(2)*1000) - params.epoch(2)*1000 -1;
%tpz = length(-1000:1500) - 1500 
baseidx = [timepoint_zero + baseline_window(1) timepoint_zero + baseline_window(2)];

%%
time = params.epoch(1):1/params.srate:params.epoch(2);
eval(sprintf('ntimes = size(block%d.time{1,1},2);', bb));
eval(sprintf('chanlabels = block%d.label;', bb));
%% work individually for for hand and mobile 
nchans = length(chanlabels);  % number of channels does not change across blocks 

for cc =1:nchans
    
    % put data in matrices
    for bb = block2keep'
        eval(sprintf('ntrials = length(block%d.trial);', bb));
        for nt = 1:ntrials
            eval(sprintf('block%d_dat(nt,:) = block%d.trial{1,nt}(cc,:) ;', bb, bb));
        end
    end

    % loop over conditions
for cond = 1:length(params.conditions);
    
    % compute number of trials per condition
    if cond == 1
        eval(sprintf('ntrials(cond) = length(block%d.trial) + length(block%d.trial);', block2keep(1), block2keep(2)));
    elseif cond == 2
        eval(sprintf('ntrials(cond) = length(block%d.trial) + length(block%d.trial);', block2keep(3), block2keep(4)));
    end
        alldata = zeros(ntrials(cond), ntimes);

    % create a matrix for each condition (each block, each channel)
    if cond == 1
        eval(sprintf('alldata(:,:) = [block%d_dat; block%d_dat] ;', block2keep(1), block2keep(2)));
    elseif cond == 2
        eval(sprintf('alldata(:,:) = [block%d_dat; block%d_dat] ;', block2keep(3), block2keep(4)));
    end
    
        
%% hilbert filter parameters
% initialize output time-frequency data
tf = zeros(1,length(time));

%reshape data matrix in concatenating all the trials in a row 
data_resh = zeros(1,ntrials(cond)*ntimes);
sizedatresh = size(data_resh);
data_resh = reshape(alldata(:,:)',[sizedatresh]);

% apply the filter kernel to the signal
filtsig = zeros(1,ntrials(cond)*ntimes);
filtsig(cc,:) = filtfilt(filtkern,1,data_resh);
% hilbert transform
hilbfiltsig = zeros(1,ntrials(cond)*ntimes);
hilbfiltsig = hilbert(filtsig(cc,:));

as_resh = zeros(ntrials(cond), ntimes);
as_resh_abs = zeros(ntrials(cond), ntimes);
as_resh = reshape(hilbfiltsig,ntimes, ntrials(cond))';
as_resh_abs(:,:) = abs(as_resh).^2;
% 
    %% power normalization
    
% power normalization within reach trial (afterwards, I can average
        % the frequency together)
        % baseline is computed as the mean baseline of the condition
        %%
    activity = zeros(ntrials(cond), ntimes);    
    tfDB = zeros(ntrials(cond), ntimes);   
    as_resh_overtr = mean(as_resh_abs,1);
    baseline_ovtr = mean(as_resh_overtr(1,baseidx(1):baseidx(2)),2); % average baselines of the trials 
    
    for iii = 1:ntrials(cond);
    activity(iii,:) = as_resh_abs(iii,:);
    baselineeachtr = mean(activity(iii,baseidx(1):baseidx(2)),2);
    baselines(iii,1) = baselineeachtr;
    % put in DB scale
    tfDB(iii,:) = 10*log10( activity(iii,:) ./baselineeachtr);
    end
    
    
%% create files to be saved later 

eval(sprintf('tfhilb_raw_%s(cc,ff,:,:) = as_resh_abs(:,:); ', params.conditions{cond}));
eval(sprintf('tfhilb_DB_%s(cc,ff,:,:) = tfDB(:,:); ',  params.conditions{cond}));


end 
%% plot power averaged over trials and frequencies 
f = figure('visible','off');
options.handle     = f;

for cond = 1:length(params.conditions);
eval(sprintf('data2plot_%s = zeros(ntrials(cond), ntimes);', params.conditions{cond}));
eval(sprintf('data2plot_%s(:, :)  = squeeze(tfhilb_DB_%s(cc,ff,:,:));',params.conditions{cond}, params.conditions{cond}));

% average over trials
eval(sprintf('data2plot_av_%s = mean(data2plot_%s,1);', params.conditions{cond}, params.conditions{cond}));
eval(sprintf('data2plot_std_%s = std(data2plot_%s,1)'';', params.conditions{cond}, params.conditions{cond}));
% smooth 
sparam = params.smoothpar(ff);
eval(sprintf('data2plot_av_smooth_%s = smooth(data2plot_av_%s, sparam, ''moving'');', params.conditions{cond}, params.conditions{cond}));
eval(sprintf('data2plot_std_smooth_%s = smooth(data2plot_std_%s, sparam, ''moving'');', params.conditions{cond}, params.conditions{cond}));
lims = [-600 1000]; % around timezero 
plot_tw = [timepoint_zero + lims(1) timepoint_zero + lims(2)];  % specify temporal window to be plotted
if ff == 1 % beta
    ylimits = [-5 3];
elseif ff ==2 % HG
    ylimits = [-5 3];
end 
colors_area = [[102 153 255]./255; [255 53 0]./255]; % red and blue
colors_line = [[51 0 255]./255; [204 0 0 ]./255];
% options for shaded area error plot
options.alpha      = 0.2;
options.line_width = 1;
options.error      = 'sem'; 
options.ntrials     = ntrials(cond);
options.color_area = colors_area(cond,:);
options.color_line = colors_line(cond,:);

eval(sprintf('plot_areaerrorbar(data2plot_av_smooth_%s(plot_tw(1):plot_tw(2))'', data2plot_std_smooth_%s(plot_tw(1):plot_tw(2))'', options );', params.conditions{cond}, params.conditions{cond}));
hold on
% plot timewindow lines
xst = [timepoint_zero - plot_tw(1) timepoint_zero - plot_tw(1) ]; % vertical line to signal time 0
xnd = [xst(1) + 600 xst(1) + 600  ]; % vertical line to signal end of time window considered in the analysis (600 ms after the stimulus onset)
y = ylimits;
plot(xst, y,'LineWidth', 1,'LineStyle', '--', 'Color', 'r', 'HandleVisibility','off');
plot(xnd, y,'LineWidth', 1,'LineStyle', '--', 'Color', 'r', 'HandleVisibility','off');
end 
% other options
xlim([0 plot_tw(2)-plot_tw(1)]);
xticklabels([lims(1):200:lims(2)]);
xlabel('msec');
[leg,labelhandles,outH,outM] =legend('Hand','Mobile');%, 'Location')%, legendpos);
legendpos = 'northeast';
title(['average ' params.targetfrex{ff} ' power ' chanlabels{cc}]);
ylabel(['normalized ' params.targetfrex{ff} ' power (dB)']);
ylim([ylimits(1) ylimits(2)]);
saveas(gca, [params.outpathhilb_plots '/' params.targetfrex{ff} '/averagepower/'  params.ssID '_Average_' params.targetfrex{ff} '_power_' chanlabels{cc}], 'png' );
close 

end  
end 

%% save data
        % matrix nchan x frequency range  x  trials x times (for each freq
        % range) 
        for cond = 1:length(params.conditions);
    eval(sprintf('save([params.OutPath ''hilbTFraw_'' params.conditions{cond}], ''tfhilb_raw_%s''); ', params.conditions{cond}));
    eval(sprintf('save([params.OutPath ''hilbTFpc_'' params.conditions{cond}], ''tfhilb_PC_%s''); ', params.conditions{cond}));
        end
% 
% 
 
