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
clusters = fieldnames(Clust);
%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('//home/andrea/Documents/MatlabFunctions/MVPA-Light/startup/')
startup_MVPA_Light

%--------------------------------------------------------------------------
% Files to work with for Odor D Night
%--------------------------------------------------------------------------

% For Odor D Night
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

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


    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [0 20];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    
    for cluster = 1:numel(clusters)


        cfg                         = [];
        cfg.channel                 = Clust.(clusters{cluster});
        Time_Freq_Odor_Clust        = ft_selectdata(cfg, Time_Freq_Odor);
        Time_Freq_Vehicle_Clust     = ft_selectdata(cfg, Time_Freq_Vehicle);

        clabel = [ones(1,size(Time_Freq_Odor_Clust.powspctrm,1)),...
        ones(1,size(Time_Freq_Vehicle_Clust.powspctrm,1))+1];

        cfg                             = [];
        cfg.classifier                  = 'svm';
        cfg.metric                      = {'auc'};
        cfg.repeat                      = 2;
        cfg.k                           = 5;
        cfg.sample_dimension            = 1;
        cfg.feature_dimension           = [2,4];
        cfg.generalization_dimension    = 3;
        cfg.dimension_names             = {'samples','channels','frequencies','time points'};

        [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
            cat(1,Time_Freq_Odor_Clust.powspctrm,...
            Time_Freq_Vehicle_Clust.powspctrm), clabel);
        
        cf_Generalization.(clusters{cluster}){subj} = cf_freqxfreq;
        resultGeneralization.(clusters{cluster}){subj}  = result_freq;
       
    end

    F = Time_Freq_Odor.freq;
    t = Time_Freq_Odor.time;
    
    savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
    save(strcat(savepath,'OdorDvsVehicle_FreqGeneralization_to20_clusters'),...
        'cf_Generalization','resultGeneralization','F','t','-v7.3');
end



%--------------------------------------------------------------------------
% Files to work with for Odor M Night
%--------------------------------------------------------------------------

% For Odor D Night
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

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


    %----------------------------------------------------------------------
    % Select Time of interest
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [0 20];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    
    for cluster = 1:numel(clusters)


        cfg                         = [];
        cfg.channel                 = Clust.(clusters{cluster});
        Time_Freq_Odor_Clust        = ft_selectdata(cfg, Time_Freq_Odor);
        Time_Freq_Vehicle_Clust     = ft_selectdata(cfg, Time_Freq_Vehicle);

        clabel = [ones(1,size(Time_Freq_Odor_Clust.powspctrm,1)),...
        ones(1,size(Time_Freq_Vehicle_Clust.powspctrm,1))+1];

        cfg                             = [];
        cfg.classifier                  = 'svm';
        cfg.metric                      = {'auc'};
        cfg.repeat                      = 2;
        cfg.k                           = 5;
        cfg.sample_dimension            = 1;
        cfg.feature_dimension           = [2,4];
        cfg.generalization_dimension    = 3;
        cfg.dimension_names             = {'samples','channels','frequencies','time points'};

        [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
            cat(1,Time_Freq_Odor_Clust.powspctrm,...
            Time_Freq_Vehicle_Clust.powspctrm), clabel);
        
        cf_Generalization.(clusters{cluster}){subj} = cf_freqxfreq;
        resultGeneralization.(clusters{cluster}){subj}  = result_freq;
       
    end

    F = Time_Freq_Odor.freq;
    t = Time_Freq_Odor.time;
    
    savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
    save(strcat(savepath,'OdorMvsVehicle_FreqGeneralization_to20_clusters'),...
        'cf_Generalization','resultGeneralization','F','t','-v7.3');
end


%--------------------------------------------------------------------------
% Files to work with for both nights
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
    
    for cluster = 1:numel(clusters)

        cfg                         = [];
        cfg.channel                 = Clust.(clusters{cluster});
        Time_Freq_OdorD_Clust       = ft_selectdata(cfg, Time_Freq_OdorD);
        Time_Freq_OdorM_Clust       = ft_selectdata(cfg, Time_Freq_OdorM);
 
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
        cfg.dimension_names             = {'samples','channels','frequencies','time points'};

        [cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
            cat(1,Time_Freq_OdorD_Clust.powspctrm,...
            Time_Freq_OdorM_Clust.powspctrm), clabel);
            
        cf_Generalization.(clusters{cluster}){subj} = cf_freqxfreq;
        
        resultGeneralization.(clusters{cluster}){subj}  = result_freq;  
    end 
    F = Time_Freq_Odor.freq;
    t = Time_Freq_Odor.time;

    savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';
    save(strcat(savepath,'OdorDvsOdorM_FreqGeneralization_to20_clusters'),...
        'cf_Generalization','resultGeneralization','F','t','-v7.3');
    
end

