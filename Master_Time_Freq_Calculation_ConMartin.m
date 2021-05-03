addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

s_tstep = 0.1; % try with 0.005
s_fstep = 0.1; % 0.005

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Organize files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------

%filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Filtered-David/CueNight/';

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    

    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Cue.label);
        Time_Freq_OdorD_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
        
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Vehicle.label);
        Time_Freq_VehicleD_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));
        
    end
     
end

save('Time_Freq_Clust_ToPlot_DNight','Time_Freq_OdorD_clust','Time_Freq_VehicleD_clust')

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

% filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

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

% save('Time_Freq_Clust_ToPlot_MNight','Time_Freq_OdorM_clust','Time_Freq_VehicleM_clust')

% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       Calculation of power change in specific bands
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Time_Freq = Time_Freq_Odor;
% 
% 
% SpindleBand     = [10 18];
% DeltaBand       = [0.5 2];
% ThetaBand       = [4 8];
% BetaBand        = [18 30];
% time = Time_Freq.time ;
% subjects = 1:23;%numel(filesOdor);
% 
% 
% SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
%     Time_Freq.freq<=SpindleBand(2));
% DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
%     Time_Freq.freq<=DeltaBand(2));
% ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
%     Time_Freq.freq<=ThetaBand(2));
% BetaIdx = find(Time_Freq.freq>=BetaBand(1) &...
%     Time_Freq.freq<=BetaBand(2));
% 
% 
% %--------------------------------------------------------------------------
% % For Declarative Associated Odor Night
% %--------------------------------------------------------------------------
% for cluster = 1:numel(clusters)
%     % -----------------For Cue----------------------------------
%     v_Spindle_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
%     v_Delta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
%     v_Theta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
%     v_Beta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
% 
%     % -----------------For Vehicle----------------------------------
%     v_Spindle_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
%     v_Delta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
%     v_Theta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
%     v_Beta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
% 
% end
% 
% %--------------------------------------------------------------------------
% % For Motor Associated Odor Night
% %--------------------------------------------------------------------------
% 
% for cluster = 1:numel(clusters)
%     % -----------------For Cue----------------------------------
%     v_Spindle_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
%     v_Delta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
%     v_Theta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
%     v_Beta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
% 
%     % -----------------For Vehicle----------------------------------
%     v_Spindle_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
%     v_Delta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
%     v_Theta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
%     v_Beta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
% 
% end
% 
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       PLOT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
% addpath('./Scripts_Wilc')
% 
% %p_plot_TimeFreq_All
% tf = [];
% 
% for cluster = 1:numel(clusters)
%     
%     figure
%     y_lims          = [];
%     x_lims_parcial  = [-5 25];
%     time_parcial    = Time_Freq_Odor.time;
%     total_subplots  = 4;
%     count           = 1;
%     frequencies     = Time_Freq_Odor.freq;
%     
%     TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
%     TF_VehicleD = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster}),1));
%     
%     
%     TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
%     TF_VehicleM = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster}),1));
%     
%     bottom  = min(min(TF_OdorD(:)),min(TF_VehicleD(:)),...
%         min(TF_OdorM(:)),min(TF_VehicleM(:)));
%     top     = max(max(TF_OdorD(:)),max(TF_VehicleD(:)),...
%         max(TF_OdorM(:)),max(TF_VehicleM(:)));
%     
%     %----------------------------------------------------------------------
%     % For Declarative Associated Odor Night
%     %----------------------------------------------------------------------
% 
%     tf(count) = subplot(total_subplots,1,count);
%     f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
%     xlim(x_lims_parcial)
%     colormap(parula)
%     caxis manual
%     caxis([bottom top]);
%     colorbar
%     title('Odor D')
%     
%  
%     count = count+1;
%     tf(count) = subplot(total_subplots,1,count);
%     f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
%     xlim(x_lims_parcial)
%     colormap(parula)
%     caxis manual
%     caxis([bottom top]);
%     colorbar
%     title('Vehicle M')
%     
%     %----------------------------------------------------------------------
%     % For Motor Associated Odor Night
%     %----------------------------------------------------------------------
% 
%     tf(count) = subplot(total_subplots,1,count);
%     f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
%     xlim(x_lims_parcial)
%     colormap(parula)
%     caxis manual
%     caxis([bottom top]);
%     colorbar
%     title('Odor M')
%     
%     count = count+1;
%     tf(count) = subplot(total_subplots,1,count);
%     f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
%     xlim(x_lims_parcial)
%     colormap(parula)
%     caxis manual
%     caxis([bottom top]);
%     colorbar
%     title('Vehicle M')
%     
%     linkaxes(tf,'x')
%     
%     %set(gcf,'position',[81,35,700,926])
%     
%     l = sgtitle(clusters{cluster});
%     set(l, 'Interpreter', 'none')
%     
%     %saveas(gcf,strcat(clusters{cluster},'Final_OdorD.png'))
%     
%     
% %     v_time = Time_Freq_Odor.time ;
% %     ValidTimeOdorD = v_time(~isnan(Time_Freq_OdorD_clust.(clusters{cluster})(1,1,:)));
% %     ValidTimeVehicleD = v_time(~isnan(Time_Freq_VehicleD_clust.(clusters{cluster})(1,1,:)));
% %     
% %     ValidTime = intersect(ValidTimeOdorD,ValidTimeVehicleD);
% %     start_sample = ceil((ValidTime(1)-v_time(1))/s_tstep)+1;
% %     end_sample = floor((ValidTime(end)-v_time(1))/s_tstep);
% 
% % start_sample = 1;
% % end_sample = length(v_time);
% 
% %     count = count+1;
% %     subplot(total_subplots,1,count)
% %     f_nonParametricTest('Beta',...
% %         ' ',' ','Odor D','Vehicle D',...
% %         v_Beta_OdorD.(clusters{cluster}),...
% %         v_Beta_VehicleD.(clusters{cluster}),...
% %         v_time,'-b',start_sample,end_sample);
% %     
% %     count = count+1;
% %     subplot(total_subplots,1,count)
% %     f_nonParametricTest('Spindle',...
% %         ' ',' ','Odor D','Vehicle D',...
% %         v_Spindle_OdorD.(clusters{cluster}),...
% %         v_Spindle_VehicleD.(clusters{cluster}),...
% %         v_time,'-b',start_sample,end_sample);
% %     
% %     count = count+1;
% %     subplot(total_subplots,1,count)
% %     f_nonParametricTest('Theta',...
% %         ' ',' ','Odor D','Vehicle D',...
% %         v_Theta_OdorD.(clusters{cluster}),...
% %         v_Theta_VehicleD.(clusters{cluster}),...
% %         v_time,'-b',start_sample,end_sample);
% %     
% %     count = count+1;
% %     tf(count) = subplot(total_subplots,1,count);
% %     f_nonParametricTest('Delta',...
% %         'Time(sec)',' ','Odor D','Vehicle D',...
% %         v_Delta_OdorD.(clusters{cluster}),...
% %         v_Delta_VehicleD.(clusters{cluster}),...
% %         v_time,'-r',start_sample,end_sample);
%     
% 
% 
% end
% 
%  