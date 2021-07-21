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

%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('/home/andrea/Documents/MatlabFunctions/MVPA-Light/startup/')
startup_MVPA_Light

% %--------------------------------------------------------------------------
% % Calculating for D Night
% %--------------------------------------------------------------------------

% filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

% filesOdor = dir(strcat(filepath,'*Odor.mat'));
% filesVehicle = dir(strcat(filepath,'*Vehicle.mat'));
% %% MVPA TF analysis
% % 
% subjects = 1:numel(filesOdor);

% disp(strcat('Sujeto: ',num2str(1)))
    
% % Load Data
% load(strcat(filepath,filesOdor(1).name));
% load(strcat(filepath,filesVehicle(1).name));
% %----------------------------------------------------------------------
% % Select Time of interest
% %----------------------------------------------------------------------
% cfg = [];
% cfg.latency = [-1 5];
% Time_Freq_OdorAll = ft_selectdata(cfg, Time_Freq_Odor);
% Time_Freq_VehicleAll = ft_selectdata(cfg, Time_Freq_Vehicle);

% for subj = 2:subjects(end)
%     % display current subject name
%     disp(strcat('Sujeto: ',num2str(subj)))
    
%     % Load Data
%     load(strcat(filepath,filesOdor(subj).name));
%     load(strcat(filepath,filesVehicle(subj).name));

%     %----------------------------------------------------------------------
%     % Select Time of interest
%     %----------------------------------------------------------------------
%     cfg = [];
%     cfg.latency = [-1 5];
%     Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Odor);
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);
    
%     cfg = [];
%     %cfg.appenddim = 1;
%     Time_Freq_OdorAll = ft_appendfreq(cfg,Time_Freq_OdorAll,Time_Freq_Odor);
%     Time_Freq_VehicleAll = ft_appendfreq(cfg,Time_Freq_VehicleAll,Time_Freq_Vehicle);
% end

% clabel = [ones(1,size(Time_Freq_OdorAll.powspctrm,1)),...
% ones(1,size(Time_Freq_VehicleAll.powspctrm,1))+1];

% cfg                             = [];
% cfg.classifier                  = 'svm';
% cfg.metric                      = {'auc'};
% cfg.repeat                      = 10;
% cfg.k                           = 5;
% cfg.sample_dimension            = 1;
% cfg.feature_dimension           = [2,4];
% cfg.generalization_dimension    = 3;
% cfg.dimension_names             = {'samples','channels','frequencies','time points'};

% [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
%     cat(1,Time_Freq_OdorAll.powspctrm,...
%     Time_Freq_VehicleAll.powspctrm), clabel);

% cf_Generalization= cf_freqxfreq;
% resultGeneralization = result_freq;  

% savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
% save(strcat(savepath,'OdorDvsVehicle_FreqGeneralization_[-1,5]_SuperSubj'),'cf_Generalization','resultGeneralization');

% %--------------------------------------------------------------------------
% % Calculating for M Night
% %--------------------------------------------------------------------------

% filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

% filesOdor = dir(strcat(filepath,'*Odor.mat'));
% filesVehicle = dir(strcat(filepath,'*Vehicle.mat'));
% %% MVPA TF analysis
% % 
% subjects = 1:numel(filesOdor);

% disp(strcat('Sujeto: ',num2str(1)))
    
% % Load Data
% load(strcat(filepath,filesOdor(1).name));
% load(strcat(filepath,filesVehicle(1).name));
% %----------------------------------------------------------------------
% % Select Time of interest
% %----------------------------------------------------------------------
% cfg = [];
% cfg.latency = [-1 5];
% Time_Freq_OdorAll = ft_selectdata(cfg, Time_Freq_Odor);
% Time_Freq_VehicleAll = ft_selectdata(cfg, Time_Freq_Vehicle);

% for subj = 2:subjects(end)
%     % display current subject name
%     disp(strcat('Sujeto: ',num2str(subj)))
    
%     % Load Data
%     load(strcat(filepath,filesOdor(subj).name));
%     load(strcat(filepath,filesVehicle(subj).name));


%     %----------------------------------------------------------------------
%     % Select Time of interest
%     %----------------------------------------------------------------------
%     cfg = [];
%     cfg.latency = [-1 5];
%     Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Odor);
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);
    
%     cfg = [];
%     %cfg.appenddim = 1;
%     Time_Freq_OdorAll = ft_appendfreq(cfg,Time_Freq_OdorAll,Time_Freq_Odor);
%     Time_Freq_VehicleAll = ft_appendfreq(cfg,Time_Freq_VehicleAll,Time_Freq_Vehicle);
% end

% clabel = [ones(1,size(Time_Freq_OdorAll.powspctrm,1)),...
% ones(1,size(Time_Freq_VehicleAll.powspctrm,1))+1];

% cfg                             = [];
% cfg.classifier                  = 'svm';
% cfg.metric                      = {'auc'};
% cfg.repeat                      = 10;
% cfg.k                           = 5;
% cfg.sample_dimension            = 1;
% cfg.feature_dimension           = [2,4];
% cfg.generalization_dimension    = 3;
% cfg.dimension_names             = {'samples','channels','frequencies','time points'};

% [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
%     cat(1,Time_Freq_OdorAll.powspctrm,...
%     Time_Freq_VehicleAll.powspctrm), clabel);

% cf_Generalization= cf_freqxfreq;
% resultGeneralization = result_freq;  

% savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
% save(strcat(savepath,'OdorMvsVehicle_FreqGeneralization_[-1,5]_SuperSubj'),'cf_Generalization','resultGeneralization');


%--------------------------------------------------------------------------
% Calculating for Both Nights
%--------------------------------------------------------------------------

filepathD = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';
filepathM = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

filesOdorD = dir(strcat(filepathD,'*Odor.mat'));
filesOdorM = dir(strcat(filepathM,'*Odor.mat'));
%% MVPA TF analysis
% 
subjects = 1:numel(filesOdorD);

disp(strcat('Sujeto: ',num2str(1)))
    
% Load Data
Time_Freq_OdorD = load(strcat(filepathD,filesOdorD(1).name));
Time_Freq_OdorM = load(strcat(filepathM,filesOdorM(1).name));
%----------------------------------------------------------------------
% Select Time of interest
%----------------------------------------------------------------------
cfg = [];
cfg.latency = [-1 5];
Time_Freq_OdorD_All = ft_selectdata(cfg, Time_Freq_OdorD.Time_Freq_Odor);
Time_Freq_OdorM_All = ft_selectdata(cfg, Time_Freq_OdorM).Time_Freq_Odor;

for subj = 2:subjects(end)
    % display current subject name
    disp(strcat('Sujeto: ',num2str(subj)))
    
    % Load Data
    Time_Freq_OdorD = load(strcat(filepathD,filesOdorD(subj).name));
    Time_Freq_OdorM = load(strcat(filepathM,filesOdorM(subj).name));


    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-1 5];
    Time_Freq_OdorD = ft_selectdata(cfg, Time_Freq_OdorD.Time_Freq_Odor);
    Time_Freq_OdorM = ft_selectdata(cfg, Time_Freq_OdorM.Time_Freq_Odor);
    
    cfg = [];
    cfg.appenddim = 'rpt';
    Time_Freq_OdorD_All = ft_appendfreq(cfg,Time_Freq_OdorD_All,Time_Freq_OdorD);
    Time_Freq_OdorM_All = ft_appendfreq(cfg,Time_Freq_OdorM_All,Time_Freq_OdorM);
end

clabel = [ones(1,size(Time_Freq_OdorD_All.powspctrm,1)),...
ones(1,size(Time_Freq_OdorM_All.powspctrm,1))+1];

cfg                             = [];
cfg.classifier                  = 'svm';
cfg.metric                      = {'auc'};
cfg.repeat                      = 10;
cfg.k                           = 5;
cfg.sample_dimension            = 1;
cfg.feature_dimension           = [2,4];
cfg.generalization_dimension    = 3;
cfg.dimension_names             = {'samples','channels','frequencies','time points'};

[cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
    cat(1,Time_Freq_OdorD_All.powspctrm,...
    Time_Freq_OdorM_All.powspctrm), clabel);

cf_Generalization= cf_freqxfreq;
resultGeneralization = result_freq;  

savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
save(strcat(savepath,'OdorDvsOdorM_FreqGeneralization_[-1,5]_SuperSubj'),'cf_Generalization','resultGeneralization');