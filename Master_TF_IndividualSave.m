% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/') %server
addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip\') %Windows

ft_defaults
% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')%server
addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip\qsub') %Windows
ft_warning off

%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.5; % try with 0.005
s_fstep = 0.5; % 0.005
% cycles  = 12;



cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  = 0.5:s_fstep:20;


% cycles = zeros(size(cfg_Tf.foi));
% cycles(1:36) = 4; %frequencies until 4Hz
% cycles(37:76) = 6; %frequencies from 4 to 8 Hz
% cycles(77:116) = 8; %frequencies from 8 to 12 Hz
% cycles(117:156) = 10; %frequencies from 12 to 16 Hz
% cycles(157:196) = 12; %frequencies from 16 to 20 Hz


cycles  = 4:8/numel(cfg_Tf.foi):12;
cfg_Tf.width                = cycles;

% v_timeWindows               = 5:(-4.5/numel(cfg_Tf.foi)):0.5;
% cfg_Tf.t_ftimwin            = v_timeWindows;

cfg_Tf.toi                  = -12:s_tstep:72; 
% toi this is extended before and after to deal with border effect of wavelet
cfg_Tf.keeptrials           = 'yes';

warning off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/DNight/';%server
% savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';%server


filepath = 'D:\Thesis_Publication\Epoched_90SecTrial_MastoidRef-Interp_NewEpochs\DNight\';%Windows
savepath = 'D:\TF_Calculation_90SecTrial_New\DNight\';%Windows

files = dir(strcat(filepath,'*.set'));


for subj = 1%:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
%     addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))%server
    addpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))%Windows
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    rmpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))%Windows
%     rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))%server

    
    addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip\external\eeglab\') %windows
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw','none');
    
    
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-5 60];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear EEGOdor dataOdor
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-5 60];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    clear Time_Freq
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [25 60];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(savepath,files(subj).name(1:12),'TF_DN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_DN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Motor Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/MNight/';%server
% savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';%server

filepath = 'D:\Thesis_Publication\Epoched_90SecTrial_MastoidRef-Interp_NewEpochs\MNight\';%Windows
savepath = 'D:\TF_Calculation_90SecTrial_New\MNight\';%Windows

files = dir(strcat(filepath,'*.set'));


for subj = 7:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
%     addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))%server
    addpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))%Windows
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    rmpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))%Windows
%     rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))%server

    
    addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip\external\eeglab\') %windows
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw','none');
    

    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    clear EEGOdor dataOdor
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-5 60];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear Time_Freq_DA_Temp
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-5 60];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [25 60];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(savepath,files(subj).name(1:12),'TF_MN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_MN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end


%% Plot Individual TF

figure
y_lims                  = [];
x_lims_parcial          = [-5 25];
time_parcial            = Time_Freq_Odor.time;
frequencies             = Time_Freq_Odor.freq;
v_xlim                  = [-5 25];
v_xticks                = [-5 0 5 10 15 20 25];

TF = mean(Time_Freq_Odor.powspctrm,1);
TF = squeeze(TF);
TF = squeeze(mean(TF,1));
max_TF = abs(TF);
max_TF = max(max_TF(:));


f_ImageMatrix(TF,time_parcial,frequencies,y_lims)
xlim(v_xlim)
ylim([0 20])
colormap(parula)
caxis manual
caxis([-max_TF max_TF]);
title('TF')
xlabel('')
hold on
plot([0,15],[0.25,0.25],'k','LineWidth',3);
hold off


