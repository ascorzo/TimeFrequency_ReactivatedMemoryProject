
ft_defaults
load('clusterChans.mat')
ft_warning off
%---------------------------------------
%Time frequency parameters
%---------------------------------------
cfg_Tf                 = [];
cfg_Tf.method          = 'mtmconvol';
cfg_Tf.taper           = 'hanning';
cfg_Tf.output          = 'pow';
cfg_Tf.pad             = 'nextpow2';
cfg_Tf.foi             =  0.5:0.25:30;
cfg_Tf.t_ftimwin       = ones(length(cfg_Tf.foi),1);   % length of time window = 0.5 sec
cfg_Tf.toi             = -15:0.25:30;
cfg_Tf.keeptrials      = 'yes';


%------------------------------------
% Parameters for Selection of data based on channels of interest 
%------------------------------------
cfg_Sel = [];
cfg_Sel.channel = [frontal_channels;central_channels];

%------------------------------------
% Parameters for Baseline Correction 
%------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-5 0];
cfg_Bas.baselinetype        = 'relative';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Cue Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = 'D:\GermanData\DATA\RawData\preProcessing\Epochs\CueNight\';

filesOdor = dir(strcat(filepath,'*Odor*.set'));
filesSham = dir(strcat(filepath,'*Vehicle*.set'));


for subj = 1:numel(filesOdor)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))

    %--------- Load Data ---------------------------------------
    addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
    EEGOdor = pop_loadset(strcat(filepath,filesOdor(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    rmpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'));
    
 
    %-----------Selection of data Channels--------------------------------- 
    dataOdor = ft_selectdata(cfg_Sel,dataOdor);
    
    %-----------Time-Frequency Calculation---------------------------------
    Time_Freq_Cue{subj}   = ft_freqanalysis(cfg_Tf, dataOdor);
    
    %-----------Baseline Correction ---------------------------------------
    Time_Freq_Cue_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq_Cue{subj});
    
    %-----------Mean of all trials, each channel --------------------------
    Time_Freq_Cue_Mean = Time_Freq_Cue_Baseline;
    Time_Freq_Cue_Mean.powspctrm = squeeze(mean(Time_Freq_Cue_Mean.powspctrm,1));
    
    %-----------Combine Mean of all frontal channels-----------------------
    [~,~,ind2] = intersect(frontal_channels,Time_Freq_Cue_Mean.label);
    Time_Freq_Cue_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Cue_Mean.powspctrm(ind2,:,:),1));
    
    %-----------Combine Mean of all central channels-----------------------
    [~,~,ind2] = intersect(central_channels,Time_Freq_Cue_Mean.label);
    Time_Freq_Cue_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Cue_Mean.powspctrm(ind2,:,:),1));
    
    %______________________________________________________________________
    %
    %Sham Odor
    %______________________________________________________________________ 

    disp(strcat('Sham'))

    %--------- Load Data --------------------------------------------------
    addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
    EEGSham = pop_loadset(strcat(filepath,filesSham(subj).name));
    dataSham = eeglab2fieldtrip(EEGSham,'raw');
    rmpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'));
     

    %-----------Selection of data Channels--------------------------------- 
    dataSham = ft_selectdata(cfg_Sel,dataSham);
    
    %-----------Time-Frequency Calculation---------------------------------
    Time_Freq_Sham_CN{subj}   = ft_freqanalysis(cfg_Tf, dataSham);
    
    %-----------Baseline Correction ---------------------------------------
    Time_Freq_Sham_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq_Sham{subj});
    
    %-----------Mean of all trials, each channel --------------------------
    Time_Freq_Sham_Mean = Time_Freq_Sham_Baseline;
    Time_Freq_Sham_Mean.powspctrm = squeeze(mean(Time_Freq_Sham_Mean.powspctrm,1));
    
    %-----------Combine Mean of all frontal channels-----------------------
    [~,~,ind2] = intersect(frontal_channels,Time_Freq_Sham_Mean.label);
    Time_Freq_Sham_CN_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Sham_Mean.powspctrm(ind2,:,:),1));
    
    %-----------Combine Mean of all central channels-----------------------
    [~,~,ind2] = intersect(central_channels,Time_Freq_Sham_Mean.label);
    Time_Freq_Sham_CN_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Sham_Mean.powspctrm(ind2,:,:),1));
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Placebo Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = 'D:\GermanData\DATA\RawData\preProcessing\Epochs\PlaceboNight\';

filesOdor = dir(strcat(filepath,'*Odor*.set'));
filesSham = dir(strcat(filepath,'*Vehicle*.set'));

for subj = 1:numel(filesOdor)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))

    %--------- Load Data ---------------------------------------
    addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
    EEGOdor = pop_loadset(strcat(filepath,filesOdor(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    rmpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'));
    
 
    %-----------Selection of data Channels--------------------------------- 
    dataOdor = ft_selectdata(cfg_Sel,dataOdor);
    
    %-----------Time-Frequency Calculation---------------------------------
    Time_Freq_Placebo{subj}   = ft_freqanalysis(cfg_Tf, dataOdor);
    
    %-----------Baseline Correction ---------------------------------------
     Time_Freq_Placebo_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq_Placebo{subj});
    
    %-----------Mean of all trials, each channel --------------------------
    Time_Freq_Placebo_Mean = Time_Freq_Placebo_Baseline;
    Time_Freq_Placebo_Mean.powspctrm = squeeze(mean(Time_Freq_Placebo_Mean.powspctrm,1));
    
    %-----------Combine Mean of all frontal channels-----------------------
    [~,~,ind2] = intersect(frontal_channels,Time_Freq_Placebo_Mean.label);
    Time_Freq_Placebo_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Placebo_Mean.powspctrm(ind2,:,:),1));
    
    %-----------Combine Mean of all central channels-----------------------
    [~,~,ind2] = intersect(central_channels,Time_Freq_Placebo_Mean.label);
    Time_Freq_Placebo_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Placebo_Mean.powspctrm(ind2,:,:),1));
    
    %______________________________________________________________________
    %
    %Sham Odor
    %______________________________________________________________________ 

    disp(strcat('Sham'))

    %--------- Load Data --------------------------------------------------
    addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
    EEGSham = pop_loadset(strcat(filepath,filesSham(subj).name));
    dataSham = eeglab2fieldtrip(EEGSham,'raw');
    rmpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'));
    

    %-----------Selection of data Channels--------------------------------- 
    dataSham = ft_selectdata(cfg_Sel,dataSham);
    
    %-----------Time-Frequency Calculation---------------------------------
    Time_Freq_Sham_PN{subj}   = ft_freqanalysis(cfg_Tf, dataSham);
    
    %-----------Baseline Correction ---------------------------------------
    Time_Freq_Sham_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq_Sham{subj});
    
    %-----------Mean of all trials, each channel --------------------------
    Time_Freq_Sham_Mean = Time_Freq_Sham_Baseline;
    Time_Freq_Sham_Mean.powspctrm = squeeze(mean(Time_Freq_Sham_Mean.powspctrm,1));
    
    %-----------Combine Mean of all frontal channels-----------------------
    [~,~,ind2] = intersect(frontal_channels,Time_Freq_Sham_Mean.label);
    Time_Freq_Sham_PN_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Sham_Mean.powspctrm(ind2,:,:),1));
    
    %-----------Combine Mean of all central channels-----------------------
    [~,~,ind2] = intersect(central_channels,Time_Freq_Sham_Mean.label);
    Time_Freq_Sham_PN_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Sham_Mean.powspctrm(ind2,:,:),1));
end

%% Calculation of power change in specific bands

SpindleBand     = [12 16];
DeltaBand       = [0.5 4];
time = Time_Freq_Placebo.time ;

%------------------------------------------------------
% For Cue Night
%------------------------------------------------------

SpindleIdx = find(Time_Freq_Cue_Mean.freq>=SpindleBand(1) &...
    Time_Freq_Cue_Mean.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq_Cue_Mean.freq>=DeltaBand(1) &...
    Time_Freq_Cue_Mean.freq<=DeltaBand(2));

v_Spindle_Cue_Frontal = squeeze(mean(Time_Freq_Cue_Frontal(:,SpindleIdx,:),2));
v_Spindle_Sham_CN_Frontal = squeeze(mean(Time_Freq_Sham_CN_Frontal(:,SpindleIdx,:),2));

v_Delta_Cue_Frontal = squeeze(mean(Time_Freq_Cue_Frontal(:,DeltaIdx,:),2));
v_Delta_Sham_CN_Frontal = squeeze(mean(Time_Freq_Sham_CN_Frontal(:,DeltaIdx,:),2));


v_Spindle_Cue_Central = squeeze(mean(Time_Freq_Cue_Central(:,SpindleIdx,:),2));
v_Spindle_Sham_CN_Central = squeeze(mean(Time_Freq_Sham_CN_Central(:,SpindleIdx,:),2));

v_Delta_Cue_Central = squeeze(mean(Time_Freq_Cue_Central(:,DeltaIdx,:),2));
v_Delta_Sham_CN_Central = squeeze(mean(Time_Freq_Sham_CN_Central(:,DeltaIdx,:),2));

%------------------------------------------------------
% For Placebo Night
%------------------------------------------------------
SpindleIdx = find(Time_Freq_Placebo_Mean.freq>=SpindleBand(1) &...
    Time_Freq_Placebo_Mean.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq_Placebo_Mean.freq>=DeltaBand(1) &...
    Time_Freq_Placebo_Mean.freq<=DeltaBand(2));

v_Spindle_Placebo_Frontal = squeeze(mean(Time_Freq_Placebo_Frontal(:,SpindleIdx,:),2));
v_Spindle_Sham_PN_Frontal = squeeze(mean(Time_Freq_Sham_PN_Frontal(:,SpindleIdx,:),2));

v_Delta_Placebo_Frontal = squeeze(mean(Time_Freq_Placebo_Frontal(:,DeltaIdx,:),2));
v_Delta_Sham_PN_Frontal = squeeze(mean(Time_Freq_Sham_PN_Frontal(:,DeltaIdx,:),2));


v_Spindle_Placebo_Central = squeeze(mean(Time_Freq_Placebo_Central(:,SpindleIdx,:),2));
v_Spindle_Sham_PN_Central = squeeze(mean(Time_Freq_Sham_PN_Central(:,SpindleIdx,:),2));

v_Delta_Placebo_Central = squeeze(mean(Time_Freq_Placebo_Central(:,DeltaIdx,:),2));
v_Delta_Sham_PN_Central = squeeze(mean(Time_Freq_Sham_PN_Central(:,DeltaIdx,:),2));

%clear Time_Freq_Cue_Mean Time_Freq_Sham_CN_Mean Time_Freq_Placebo_Mean Time_Freq_Sham_PN_Mean
%% plot CUE NIGHT
%----plot frontal channels----

% SPINDLE plot

figure

subplot(2,1,1)

SEM_temp_Cue = std(v_Spindle_Cue_Frontal,[],1)/sqrt(size(v_Spindle_Cue_Frontal,1));
MEAN_temp_Cue = mean(v_Spindle_Cue_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Cue-SEM_temp_Cue/2, MEAN_temp_Cue+SEM_temp_Cue/2,...
    'r','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Cue,'r');
hold on
SEM_temp_Sham_CN = std(v_Spindle_Sham_CN_Frontal,[],1)/sqrt(size(v_Spindle_Sham_CN_Frontal,1));
MEAN_temp_Sham_CN = mean(v_Spindle_Sham_CN_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Sham_CN-SEM_temp_Sham_CN/2, MEAN_temp_Sham_CN+SEM_temp_Sham_CN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_CN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (12-16Hz)')
title('Power change for frontal channels')

legend([h1 h2],{'cue' 'sham'});

% DELTA plot
subplot(2,1,2)

SEM_temp_Cue = std(v_Delta_Cue_Frontal,[],1)/sqrt(size(v_Delta_Cue_Frontal,1));
MEAN_temp_Cue = mean(v_Delta_Cue_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Cue-SEM_temp_Cue/2, MEAN_temp_Cue+SEM_temp_Cue/2,...
    'r','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Cue,'r');
hold on
SEM_temp_Sham_CN = std(v_Delta_Sham_CN_Frontal,[],1)/sqrt(size(v_Delta_Sham_CN_Frontal,1));
MEAN_temp_Sham_CN = mean(v_Delta_Sham_CN_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Sham_CN-SEM_temp_Sham_CN/2, MEAN_temp_Sham_CN+SEM_temp_Sham_CN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_CN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (0.5-4Hz)')

legend([h1 h2],{'cue' 'sham'});

%----plot central channels----

% SPINDLE plot

figure

subplot(2,1,1)

SEM_temp_Cue = std(v_Spindle_Cue_Central,[],1)/sqrt(size(v_Spindle_Cue_Central,1));
MEAN_temp_Cue = mean(v_Spindle_Cue_Central,1);
shadedplot(time, ...
    MEAN_temp_Cue-SEM_temp_Cue/2, MEAN_temp_Cue+SEM_temp_Cue/2,...
    'r','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Cue,'r');
hold on
SEM_temp_Sham_CN = std(v_Spindle_Sham_CN_Central,[],1)/sqrt(size(v_Spindle_Sham_CN_Central,1));
MEAN_temp_Sham_CN = mean(v_Spindle_Sham_CN_Central,1);
shadedplot(time, ...
    MEAN_temp_Sham_CN-SEM_temp_Sham_CN/2, MEAN_temp_Sham_CN+SEM_temp_Sham_CN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_CN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (12-16Hz)')
title('Power change for central channels')

legend([h1 h2],{'cue' 'sham'});


% DELTA plot
subplot(2,1,2)

SEM_temp_Cue = std(v_Delta_Cue_Central,[],1)/sqrt(size(v_Delta_Cue_Central,1));
MEAN_temp_Cue = mean(v_Delta_Cue_Central,1);
shadedplot(time, ...
    MEAN_temp_Cue-SEM_temp_Cue/2, MEAN_temp_Cue+SEM_temp_Cue/2,...
    'r','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Cue,'r');
hold on
SEM_temp_Sham_CN = std(v_Delta_Sham_CN_Central,[],1)/sqrt(size(v_Delta_Sham_CN_Central,1));
MEAN_temp_Sham_CN = mean(v_Delta_Sham_CN_Central,1);
shadedplot(time, ...
    MEAN_temp_Sham_CN-SEM_temp_Sham_CN/2, MEAN_temp_Sham_CN+SEM_temp_Sham_CN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_CN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (0.5-4Hz)')
legend([h1 h2],{'cue' 'sham'});

%% plot PLACEBO NIGHT
%----plot frontal channels----

% SPINDLE plot

figure

subplot(2,1,1)

SEM_temp_Placebo = std(v_Spindle_Placebo_Frontal,[],1)/sqrt(size(v_Spindle_Placebo_Frontal,1));
MEAN_temp_Placebo = mean(v_Spindle_Placebo_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Placebo-SEM_temp_Placebo/2, MEAN_temp_Placebo+SEM_temp_Placebo/2,...
    'b','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Placebo,'b');
hold on
SEM_temp_Sham_PN = std(v_Spindle_Sham_PN_Frontal,[],1)/sqrt(size(v_Spindle_Sham_PN_Frontal,1));
MEAN_temp_Sham_PN = mean(v_Spindle_Sham_PN_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Sham_PN-SEM_temp_Sham_PN/2, MEAN_temp_Sham_PN+SEM_temp_Sham_PN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_PN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (12-16Hz)')
title('Power change for frontal channels')

legend([h1 h2],{'placebo' 'sham'});

% DELTA plot
subplot(2,1,2)

SEM_temp_Placebo = std(v_Delta_Placebo_Frontal,[],1)/sqrt(size(v_Delta_Placebo_Frontal,1));
MEAN_temp_Placebo = mean(v_Delta_Placebo_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Placebo-SEM_temp_Placebo/2, MEAN_temp_Placebo+SEM_temp_Placebo/2,...
    'b','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Placebo,'b');
hold on
SEM_temp_Sham_PN = std(v_Delta_Sham_PN_Frontal,[],1)/sqrt(size(v_Delta_Sham_PN_Frontal,1));
MEAN_temp_Sham_PN = mean(v_Delta_Sham_PN_Frontal,1);
shadedplot(time, ...
    MEAN_temp_Sham_PN-SEM_temp_Sham_PN/2, MEAN_temp_Sham_PN+SEM_temp_Sham_PN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_PN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (0.5-4Hz)')
xlabel('time (s)')

legend([h1 h2],{'placebo' 'sham'});


%----plot central channels----

% SPINDLE plot

figure

subplot(2,1,1)

SEM_temp_Placebo = std(v_Spindle_Placebo_Central,[],1)/sqrt(size(v_Spindle_Placebo_Central,1));
MEAN_temp_Placebo = mean(v_Spindle_Placebo_Central,1);
shadedplot(time, ...
    MEAN_temp_Placebo-SEM_temp_Placebo/2, MEAN_temp_Placebo+SEM_temp_Placebo/2,...
    'b','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Placebo,'b');
hold on
SEM_temp_Sham_PN = std(v_Spindle_Sham_PN_Central,[],1)/sqrt(size(v_Spindle_Sham_PN_Central,1));
MEAN_temp_Sham_PN = mean(v_Spindle_Sham_PN_Central,1);
shadedplot(time, ...
    MEAN_temp_Sham_PN-SEM_temp_Sham_PN/2, MEAN_temp_Sham_PN+SEM_temp_Sham_PN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_PN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (12-16Hz)')
title('Power change for central channels')

legend([h1 h2],{'placebo' 'sham'});


% DELTA plot
subplot(2,1,2)

SEM_temp_Placebo = std(v_Delta_Placebo_Central,[],1)/sqrt(size(v_Delta_Placebo_Central,1));
MEAN_temp_Placebo = mean(v_Delta_Placebo_Central,1);
shadedplot(time, ...
    MEAN_temp_Placebo-SEM_temp_Placebo/2, MEAN_temp_Placebo+SEM_temp_Placebo/2,...
    'b','none');
alpha(.2)
hold on
h1 = plot(time,MEAN_temp_Placebo,'b');
hold on
SEM_temp_Sham_PN = std(v_Delta_Sham_PN_Central,[],1)/sqrt(size(v_Delta_Sham_PN_Central,1));
MEAN_temp_Sham_PN = mean(v_Delta_Sham_PN_Central,1);
shadedplot(time, ...
    MEAN_temp_Sham_PN-SEM_temp_Sham_PN/2, MEAN_temp_Sham_PN+SEM_temp_Sham_PN/2,...
    'k','none');
alpha(.2)
hold on
h2 = plot(time,MEAN_temp_Sham_PN,'k');
xlim([-5 30])
ylim([0.9 1.8])
ylabel('Pow rel change (0.5-4Hz)')
xlabel('time (s)')

legend([h1 h2],{'placebo' 'sham'});

