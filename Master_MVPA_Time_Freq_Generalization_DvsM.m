%% Analyze Time Frequency Calculation
% This notebook takes the Time-frequency calculation for the Reactivated connectivity 
% project, loading each file individually and calculating time series and statistics
%% Importing Fieldtrip libraries
% In this section we add the the Fieldtrip folder in the path and run ft_defaults 
% to have access to all functions of fieldtrip

addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('//home/andrea/Documents/MatlabFunctions/MVPA-Light/startup/')
startup_MVPA_Light

%--------------------------------------------------------------------------
% Files to work with
%--------------------------------------------------------------------------

filepath_DNight = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';
filepath_MNight = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

files_DNight = dir(strcat(filepath_DNight,'*Odor.mat'));
files_MNight = dir(strcat(filepath_MNight,'*Odor.mat'));
%% MVPA TF analysis
% 
subjects = 1:numel(files_DNight);

for subj = subjects
    % display current subject name
    disp(strcat('Sujeto: ',num2str(subj)))
    
    % Load Data
    Time_freq_OdorD = load(strcat(filepath_DNight,files_DNight(subj).name));
    Time_freq_OdorM = load(strcat(filepath_MNight,files_MNight(subj).name));

    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [0 20];
    Time_Freq_OdorD = ft_selectdata(cfg, Time_freq_OdorD.Time_Freq_Odor);
    Time_Freq_OdorM = ft_selectdata(cfg, Time_freq_OdorM.Time_Freq_Odor);
    
    
 
    clabel = [ones(1,size(Time_Freq_OdorD.powspctrm,1)),...
    ones(1,size(Time_Freq_OdorM.powspctrm,1))+1];

    cfg                             = [];
    cfg.classifier                  = 'svm';
    cfg.metric                      = {'auc'};
    cfg.repeat                      = 2;
    cfg.k                           = 5;
    cfg.sample_dimension            = 1;
    cfg.feature_dimension           = [2,4];
    cfg.generalization_dimension    = 3;
    cfg.dimension_names     = {'samples','channels','frequencies','time points'};

    [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
        cat(1,Time_Freq_OdorD.powspctrm,...
        Time_Freq_OdorM.powspctrm), clabel);
        
    cf_Generalization{subj} = cf_freqxfreq;
    
    resultGeneralization{subj}  = result_freq;   
    
end


for subj = 1:numel(cf_Generalization)
    
    figure
    F = Time_Freq_Odor.freq;
    mv_plot_2D(cf_Generalization{subj}, 'x', F, 'y', F)
    xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
    title('Frequency generalization using channels-x-times as features')

end
% 
result_average_Generalization = mv_combine_results(resultGeneralization, 'average');

figure
F = Time_Freq_Odor.freq;
mv_plot_2D(result_average_Generalization.perf{1}, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Frequency generalization using channels x time as features')




