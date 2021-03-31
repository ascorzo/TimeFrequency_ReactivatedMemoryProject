addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/')

ft_defaults
addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/qsub')

ft_warning off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Filtered-David/CueNight/';

filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
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
        Time_Freq_Cue_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
        
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Vehicle.label);
        Time_Freq_Vehicle_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));
        
    end
     
end



%% Calculation of power change in specific bands
Time_Freq = Time_Freq_Odor;


SpindleBand     = [10 18];
DeltaBand       = [0.5 2];
ThetaBand       = [4 8];
BetaBand        = [18 30];
time = Time_Freq.time ;
subjects = 1:numel(filesOdor);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
BetaIdx = find(Time_Freq.freq>=BetaBand(1) &...
    Time_Freq.freq<=BetaBand(2));

for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

end
%% plot CUE Night

addpath(genpath('/gpfs01/born/group/Andrea/eeglab2019_1/'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];

for cluster = 1:numel(clusters)
    
    y_lims = [];%[-0.1 0.1];
    x_lims_parcial = [-5 25];
    time_parcial = Time_Freq_Cue.time;
    %z_lims = [min(Time_Freq_Cue_Frontal(:)),max(Time_Freq_Cue_Frontal(:))];
    
    figure
    total_subplots = 3;
    count = 1;
    frequencies = Time_Freq_Cue.freq;
    
    tf(count) = subplot(total_subplots,1,count);
    f_ImageMatrix(squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(parula)
    colorbar('northoutside')
    title('TF Odor')
    
 
    count = count+1;
    tf(count) = subplot(total_subplots,1,count);
    f_ImageMatrix(squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(parula)
    colorbar('northoutside')
    title('TF Vehicle')
    
    
    v_time = Time_Freq_Cue.time ;
    ValidTimeCue = v_time(~isnan(Time_Freq_Cue_clust.(clusters{cluster})(1,1,:)));
    ValidTimeVehicle = v_time(~isnan(Time_Freq_Vehicle_clust.(clusters{cluster})(1,1,:)));
    
    ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
    start_sample = ceil((ValidTime(1)-v_time(1))/s_tstep)+1;
    end_sample = floor((ValidTime(end)-v_time(1))/s_tstep);

% start_sample = 1;
% end_sample = length(v_time);

%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Beta',...
%         ' ',' ','Odor M','Vehicle',...
%         v_Beta_Cue.(clusters{cluster}),...
%         v_Beta_Vehicle.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
%     
%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Spindle',...
%         ' ',' ','Odor M','Vehicle',...
%         v_Spindle_Cue.(clusters{cluster}),...
%         v_Spindle_Vehicle.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
%     
%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Theta',...
%         ' ',' ','Odor M','Vehicle',...
%         v_Theta_Cue.(clusters{cluster}),...
%         v_Theta_Vehicle.(clusters{cluster}),...
%         v_time,'-b',start_sample,end_sample);
    
    count = count+1;
    tf(count) = subplot(total_subplots,1,count);
    f_nonParametricTest('Delta',...
        'Time(sec)',' ','Odor D','Vehicle',...
        v_Delta_Cue.(clusters{cluster}),...
        v_Delta_Vehicle.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample);
    
    linkaxes(tf,'x')
    
    set(gcf,'position',[81,35,700,926])
    
    l = suptitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    saveas(gcf,strcat(clusters{cluster},'Final_OdorD.png'))

end

 