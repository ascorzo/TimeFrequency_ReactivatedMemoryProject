addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

s_tstep = 0.05; % try with 0.005
s_fstep = 0.05; % 0.005

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Organize files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

filename     = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/dummyfile.set';
sensors = ft_read_sens(filename);

for subj = 1:numel(filesOdor)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
    
    
    

    %-----------Mean of all trials, each channel --------------------------
%     cfg = [];
%     cfg.avgoverrpt = 'yes';
%     %cfg.latency = [-5 25];
%     Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
%     Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
%     
% 
%     Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

     
end

% save('Time_Freq_Clust_ToPlot_DNight_7Sec','Time_Freq_OdorD_clust','Time_Freq_VehicleD_clust')

% %--------------------------------------------------------------------------
% % For Motor Associated Odor Night
% %--------------------------------------------------------------------------

% filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/MNight/';

% filesOdor = dir(strcat(filepath,'*_Odor.mat'));
% filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
% p_clustersOfInterest
% clusters = fieldnames(Clust);

% for subj = 1:numel(filesOdor)

%     disp(strcat('Sujeto: ',num2str(subj)))


%     load(strcat(filepath,filesOdor(subj).name));
%     load(strcat(filepath,filesVehicle(subj).name));

%     %p_plotSubject_ConMartin
%     %-----------Mean of all trials, each channel --------------------------
%     cfg = [];
%     cfg.avgoverrpt = 'yes';
%     %cfg.latency = [-5 25];
%     Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
%     Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    

%     Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

%     %----------------------------------------------------------------------
%     %           Combine Mean of channels by clusters
%     %----------------------------------------------------------------------
    
%     for cluster = 1:numel(clusters)
%         [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Cue.label);
%         Time_Freq_OdorM_clust.(clusters{cluster})(subj,:,:) = ...
%             squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
        
%         [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Vehicle.label);
%         Time_Freq_VehicleM_clust.(clusters{cluster})(subj,:,:) = ...
%             squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));
        
%     end
     
% end

% save('Time_Freq_Clust_ToPlot_MNight_7Sec','Time_Freq_OdorM_clust','Time_Freq_VehicleM_clust')

