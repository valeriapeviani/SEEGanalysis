clear all
close all
%% this add paths of datafiles and fieldtrip toolbox
addpath('XXX')
addpath('XXX')
addpath('XXX')

%% preprocessing patient 1
params = [];
params.ssID = 'pt1';
params.nblocks = 4; % number of EDF files (corresponding to experimental blocks or conditions)
params.pathSEEG    = 'XXX\'; % path of edf file for that patient
params.outpath_plots = 'XXX\'; % where to save plots 
params.filePSY     = '';  % name of behavioural file 1 (needed to check triggers)
params.measures    = ''  % name of behaviorual file 2 (needed in Populaetrialinfo)
params.OutPath     = 'XXX\' % where to save TFA results
params.trigchan    = {'DC DC09'};
params.nonphyschan = {'-MISC E', '-DC DC01', '-DC DC02', '-DC DC03', '-DC DC04', '-DC DC09', '-MISC SpO2', '-EEG BP1', '-EEG BP2', '-EEG BP3', '-EEG BP4', '-STIM Events'}
params.unknownchan = {'-EEG cz', '-EEG pz', '-EEG c3', '-EEG p3', '-EEG c4', '-EEG p4', '-EEG X142', '-EEG X85'}
params.epoch = [-1 1.5];
params.notch = [50 100 150];
params.bandpass = [3 160];
params.targetfrex = {'beta', 'HG'};   % frequency ranges of interest
params.frex_beta = [12.5 15 25]; % minfreq, number of freq, max freq
params.frex_HG = [70 80 150]; % minfreq, number of freq, max freq
params.range_cycles =[4 10];  % wavelet number of cycles
params.srate = 1000;   % sampling rate
params.cutidx = 50; % ms to remove from the end of the trial (edge artifact)
params.smoothpar = [120 120 ]; % smoothing parameters for beta and high-gamma
params.conditions = {'cond1', 'cond2'};
params.WMchan = {'-B''7','-B''8', '-B''9','-B''12', '-C''6', '-C''7', '-C''8', '-C''9', '-C''10', '-C''11', '-C''12', '-C''13', '-C''15', '-D''3', '-D''7', '-D''8', '-D''10', '-D''11', '-D''12', '-D''13', '-E''7', '-E''8', '-E''11', '-F''2', '-F''3', '-F''4', '-F''6', '-F''7', '-F''8', ...
    '-F''9', '-F''10', '-F''11', '-F''12', '-F''13', '-I''2', '-I''3', '-I''4', '-I''5', '-I''6', '-I''7','-K''2', '-K''3', '-K''4', '-K''7', '-K''8', '-K''11','N''3', '-N''4', '-N''5', '-N''6', '-N''7', '-N''8', '-N''9', '-N''10', '-N''11', '-N''12', '-N''13','-N''14','-N''15', '-O''3', '-O''4',...
    '-O''5', '-O''6', '-P''3', '-P''4', '-P''5', '-P''6', '-P''7','-P''9', '-R''4', '-R''5','-R''10', '-S''1', '-S''6', '-S''7', '-S''8', '-S''9', '-T''1', '-T''4', '-T''5','-T''6','-U''3', '-U''5', '-W''1', '-W''4', '-W''5', '-W''7', '-W''9' '-X''1','-X''5', '-X''6', '-X''7','-X''10', '-Y''2', '-Y''3', '-Y''4', '-Y''5', '-Y''6', ...
    '-Y''7', '-Y''8', '-Y''9', '-Y''10', '-Y''11', '-Y''12', '-Y''13', '-Y''14', '-Y''15', '-Y''16', '-Y''17'}; 
params.cleanchannnames = 'yes';   % need to clean channel names?
params.pathbehav = 'XXX\' % path for behavioural data



