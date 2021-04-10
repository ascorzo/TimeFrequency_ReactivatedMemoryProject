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

%% Initial parameters
%addpath('/mnt/disk1/andrea/German_Study/');
p_clustersOfInterest


cfg_Sel = [];
cfg_Sel.channel = Clust.occipital;

%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('//home/andrea/Documents/MatlabFunctions/MVPA-Light/startup/')
startup_MVPA_Light

%--------------------------------------------------------------------------
% Files to work with
%--------------------------------------------------------------------------

% For Odor D Night
filepath_DNight = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

filesOdor_DNight = dir(strcat(filepath_DNight,'*Odor.mat'));
filesVehicle_DNight = dir(strcat(filepath_DNight,'*Vehicle.mat'));

% For Odor M Night
filepath_MNight = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

filesOdor_MNight = dir(strcat(filepath_MNight,'*Odor.mat'));
filesVehicle_MNight = dir(strcat(filepath_MNight,'*Vehicle.mat'));

%% MVPA TF analysis
% 
subjects = 1:numel(filesOdor_DNight);


for subj = subjects
    % display current subject name
    disp(strcat('Sujeto: ',num2str(subj)))
    
    % Load Data
    D_Night_Odor = load(strcat(filepath_DNight,filesOdor_DNight(subj).name));
    D_Night_Vehicle = load(strcat(filepath_DNight,filesVehicle_DNight(subj).name));

    M_Night_Odor = load(strcat(filepath_MNight,filesOdor_MNight(subj).name));
    M_Night_Vehicle = load(strcat(filepath_MNight,filesVehicle_MNight(subj).name));


    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [0 20];
    D_Night_Time_Freq_Odor = ft_selectdata(cfg, D_Night_Odor.Time_Freq_Odor);
    D_Night_Time_Freq_Vehicle = ft_selectdata(cfg, D_Night_Vehicle.Time_Freq_Vehicle);

    M_Night_Time_Freq_Odor = ft_selectdata(cfg, M_Night_Odor.Time_Freq_Odor);
    M_Night_Time_Freq_Vehicle = ft_selectdata(cfg, M_Night_Vehicle.Time_Freq_Vehicle);
    
    %----------------------------------------------------------------------
    % Balance data
    %----------------------------------------------------------------------

    s_minTrials = min([size(D_Night_Time_Freq_Odor.powspctrm,1),...
    size(D_Night_Time_Freq_Vehicle.powspctrm,1),...
    size(M_Night_Time_Freq_Odor.powspctrm,1),...
    size(M_Night_Time_Freq_Vehicle.powspctrm,1)]);

 
    clabel = [ones(1,s_minTrials),...
    ones(1,s_minTrials)+1,...
    ones(1,s_minTrials)+2,...
    ones(1,s_minTrials)+3];

    p_MVPA_TFclassifier_4Conditions
    
    cf_Generalization{subj} = cf_freqxfreq;
    
    resultGeneralization{subj}  = result_freq;

end

savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
save(strcat(savepath,'FreqGeneralization_4Conditions'),'cf_Generalization','resultGeneralization')


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
result_average_Generalization = mv_combine_results(resultGeneralization, 'average');

figure
F = Time_Freq_Odor.freq;
mv_plot_2D(result_average_Generalization.perf{1}, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Frequency generalization using channels x time as features')




