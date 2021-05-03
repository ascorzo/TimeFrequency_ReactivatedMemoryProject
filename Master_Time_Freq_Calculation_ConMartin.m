addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

s_tstep = 0.1; % try with 0.005
s_fstep = 0.15; % 0.005

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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculation of power change in specific bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Freq = Time_Freq_Odor;
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

SpindleBand     = [10 18];
DeltaBand       = [0.5 2];
Around_2_5      = [2 8];
ThetaBand       = [4 8];
BetaBand        = [18 30];
time            = Time_Freq.time;
subjects = 1:23;%numel(filesOdor);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
BetaIdx = find(Time_Freq.freq>=BetaBand(1) &...
    Time_Freq.freq<=BetaBand(2));
Around_2_5Idx = find(Time_Freq.freq>=Around_2_5(1) &...
    Time_Freq.freq<=Around_2_5(2));


%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------
for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    v_Around_2_5_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    v_Around_2_5_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

end

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    v_Around_2_5_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    v_Around_2_5_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

for cluster = 1:numel(clusters)
    
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = Time_Freq_Odor.time;
    total_subplots_row      = 3;
    total_subplots_column   = 2;
    count                   = 1;
    frequencies             = Time_Freq_Odor.freq;
    v_xlim                  = [-5 25];
    
    TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
    TF_VehicleD = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster}),1));
    
    
    TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
    TF_VehicleM = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD(:)),min(TF_VehicleD(:)),...
        min(TF_OdorM(:)),min(TF_VehicleM(:))])*0.6;
    top_TF     = max([max(TF_OdorD(:)),max(TF_VehicleD(:)),...
        max(TF_OdorM(:)),max(TF_VehicleM(:))])*0.6;
    
    bottom_Delta  = min([min(v_Delta_OdorD.(clusters{cluster})(:)),...
        min(v_Delta_VehicleD.(clusters{cluster})(:)),...
        min(v_Delta_OdorM.(clusters{cluster})(:)),...
        min(v_Delta_VehicleM.(clusters{cluster})(:))])*0.2;
    
    top_Delta     = max([max(v_Delta_OdorD.(clusters{cluster})(:)),...
        max(v_Delta_VehicleD.(clusters{cluster})(:)),...
        max(v_Delta_OdorM.(clusters{cluster})(:)),...
        max(v_Delta_VehicleM.(clusters{cluster})(:))])*0.2;
    
    start_sample = 1;
    end_sample = length(time_parcial);
    

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
    
    % -- D Night ----
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
        'Position',[0.13,0.69,0.366055226824458,0.209574269707028],...
        'PlotBoxAspectRatio',[1,0.396226399912024,0.396226399912024]);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)','Position',[-6.977238662891239,-3.349991562537163,1.000000000000014])
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor D')
    
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
         'Position',[0.53,0.69,0.359635108481264,0.209574269707028]);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor M')
    
 
    %----------------------------------------------------------------------
    % Plot TF Vehicle
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
         'Position',[0.13,0.4,0.368027613412229,0.215735294117647]);
    f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Vehicle D')
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
         'Position',[0.53,0.4,0.357662721893491,0.215735294117647]);
    f_ImageMatrix(TF_VehicleM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Vehicle M')
    
     colorbar('Position',...
         [0.92915844819818,0.399715504978663,0.013642340756455,0.500711237553343]);
     
    
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
        'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
    f_nonParametricTest('Delta',...
        'Time(sec)',' ','Odor D','Vehicle D',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_VehicleD.(clusters{cluster}),...
        time_parcial,'-r','-k',start_sample,end_sample);
    legend('Position',...
        [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])
    
    c1 =  v_Delta_OdorD.(clusters{cluster});
    c2 =  v_Delta_VehicleD.(clusters{cluster});
    
    yaxis = [nanmean(nanmean([c1;c2]))-7*nanstd(nanstd([c1;c2])), ...
        nanmean(nanmean([c1;c2]))+7*nanstd(nanstd([c1;c2]))]; % define yaxis limits

    ylim(yaxis);
 
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
        'Position',[0.53,0.12,0.357662721893491,0.215735294117647]);
    f_nonParametricTest('Delta',...
        'Time(sec)',' ','Odor M','Vehicle M',...
        v_Delta_OdorM.(clusters{cluster}),...
        v_Delta_VehicleM.(clusters{cluster}),...
        time_parcial,'-b','-k',start_sample,end_sample);
    legend('Position',...
        [0.900751386638303,0.19006499512388,0.094674554890428,0.0462304397802772]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])
    ylim(yaxis);
    
    
    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    set(gcf,'position',[190,275,1014,703])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    saveas(gcf,strcat(clusters{cluster},'.png'))

end

%     v_time = Time_Freq_Odor.time ;
%     ValidTimeOdorD = v_time(~isnan(Time_Freq_OdorD_clust.(clusters{cluster})(1,1,:)));
%     ValidTimeVehicleD = v_time(~isnan(Time_Freq_VehicleD_clust.(clusters{cluster})(1,1,:)));
%     
%     ValidTime = intersect(ValidTimeOdorD,ValidTimeVehicleD);
%     start_sample = ceil((ValidTime(1)-v_time(1))/s_tstep)+1;
%     end_sample = floor((ValidTime(end)-v_time(1))/s_tstep);

%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Beta',...
%         ' ',' ','Odor D','Vehicle D',...
%         v_Beta_OdorD.(clusters{cluster}),...
%         v_Beta_VehicleD.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
%     
%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Spindle',...
%         ' ',' ','Odor D','Vehicle D',...
%         v_Spindle_OdorD.(clusters{cluster}),...
%         v_Spindle_VehicleD.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
%     
%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Theta',...
%         ' ',' ','Odor D','Vehicle D',...
%         v_Theta_OdorD.(clusters{cluster}),...
%         v_Theta_VehicleD.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
%     

