%% Analyze Time Frequency Calculation
% This notebook takes the Time-frequency calculation for the Reactivated connectivity 
% project, loading each file individually and calculating time series and statistics
%% Importing Fieldtrip libraries
% In this section we add the the Fieldtrip folder in the path and run ft_defaults 
% to have access to all functions of fieldtrip

addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/')

ft_defaults
addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/qsub')

ft_warning off

%% Initial parameters
%addpath('/mnt/disk1/andrea/German_Study/');
p_clustersOfInterest


cfg_Sel = [];
cfg_Sel.channel = Clust.occipital;

%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('/gpfs01/born/group/Andrea/Github/MVPA-Light/startup/')
startup_MVPA_Light

%--------------------------------------------------------------------------
% Files to work with
%--------------------------------------------------------------------------

% For Odor D Night
filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorM_Night/';

% For Odor M Night
%filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_OdorD%_Night/';

filesOdor = dir(strcat(filepath,'*Odor.mat'));
filesVehicle = dir(strcat(filepath,'*Vehicle.mat'));
%% MVPA TF analysis
% 
subjects = 1:numel(filesOdor);


for subj = subjects
    % display current subject name
    disp(strcat('Sujeto: ',num2str(subj)))
    
    % Load Data
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));

%     %----------------------------------------------------------------------
%     % baseline correction
%     %----------------------------------------------------------------------
%     cfg_Bas                     = [];
%     cfg_Bas.baseline            = [-15 45];
%     cfg_Bas.baselinetype        = 'zscore';
%     Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
%     %----------------------------------------------------------------------
%     % Separate conditions
%     %----------------------------------------------------------------------
%     cfg = [];
%     cfg.latency = [-15 15];
%     Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
%     
%     cfg = [];
%     cfg.latency = [15 45];
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
%     Time_Freq_Vehicle.time = Time_Freq_Odor.time;


    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [0 20];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    
 
    clabel = [ones(1,size(Time_Freq_Odor.powspctrm,1)),...
    ones(1,size(Time_Freq_Vehicle.powspctrm,1))+1];

    p_MVPA_TFclassifier2
    
    cf_chanGeneralization{subj} = cf_freqxfreq;
    
    resultchanGeneralization{subj}  = result_freq;
    
    
%     figure
%     F = Time_Freq_Odor.freq;
%     mv_plot_2D(cf_freqxfreq, 'x', F, 'y', F)
%     xlabel('Test frequency [Hz]'), ylabel('Train frequency [Hz]')
%     title('Frequency generalization using channels-x-times as features')
%     
    
end


% for subj = 1:numel(cf_freqGeneralization)
%     
%     figure
%     F = Time_Freq_Odor.freq;
%     mv_plot_2D(cf_freqGeneralization{subj}, 'x', F, 'y', F)
%     xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
%     title('Frequency generalization using channels-x-times as features')
% 
%  
%    
% end
% % 
result_average_ChanGeneralization = mv_combine_results(resultchanGeneralization, 'average');

figure
F = 1:numel(Time_Freq_Odor.elec.label);
mv_plot_2D(result_average_ChanGeneralization.perf{1}, 'x', F, 'y', F)
xlabel('Test Channel'), ylabel('Train Channel')
title('Channel generalization using Time-x-frequencies as features')

%%
for subj = 1:numel(cf_freqGeneralization)
    for f = 1:size(cf_freqGeneralization{subj},1)
        FreqDiagonal(subj,f) = cf_freqGeneralization{subj}(f,f);
    end
    
    figure
    F =Time_Freq_Odor.freq;
    plot(F,FreqDiagonal(subj,:))
end



%%

for subj = 1:numel(cf_freqGeneralization)
    for f = 1:size(cf_freqGeneralization{subj},1)
        FreqDiagonal(subj,f) = cf_freqGeneralization{subj}(f,f);
    end
    
    figure
    F =Time_Freq_Odor.freq;
    plot(F,FreqDiagonal(subj,:))
end




