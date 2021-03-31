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
    
    [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Odor.label);
    Time_Freq_Cue_subj_clust.(clusters{cluster}) = ...
        squeeze(mean(Time_Freq_Odor.powspctrm(:,ind2,:,:),1));
    
    [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Vehicle.label);
    Time_Freq_Vehicle_subj_clust.(clusters{cluster}) = ...
        squeeze(mean(Time_Freq_Vehicle.powspctrm(:,ind2,:,:),1));
        
    % -----------------For Cue----------------------------------
    v_Spindle_Cue_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_subj_clust.(clusters{cluster})(:,SpindleIdx,:),2));
    v_Delta_Cue_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_subj_clust.(clusters{cluster})(:,DeltaIdx,:),2));
    v_Theta_Cue_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_subj_clust.(clusters{cluster})(:,ThetaIdx,:),2));
    v_Beta_Cue_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_subj_clust.(clusters{cluster})(:,BetaIdx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_Vehicle_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_subj_clust.(clusters{cluster})(:,SpindleIdx,:),2));
    v_Delta_Vehicle_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_subj_clust.(clusters{cluster})(:,DeltaIdx,:),2));
    v_Theta_Vehicle_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_subj_clust.(clusters{cluster})(:,ThetaIdx,:),2));
    v_Beta_Vehicle_subj.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_subj_clust.(clusters{cluster})(:,BetaIdx,:),2));

end
%% plot CUE Night

addpath(genpath('/gpfs01/born/group/Andrea/eeglab2019_1/'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All

for cluster = 1:numel(clusters)
    
    y_lims = [];%[-0.1 0.1];
    x_lims_parcial = [-5 30];
    time_parcial = Time_Freq.time;
    %z_lims = [min(Time_Freq_Cue_Frontal(:)),max(Time_Freq_Cue_Frontal(:))];
    
    figure
    total_subplots = 6;
    count = 1;
    frequencies = Time_Freq.freq;
    
    subplot(total_subplots,1,count)
    f_ImageMatrix(squeeze(mean(Time_Freq_Cue_subj_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(pink)
    %colorbar
    title('TF Odor')
 
    count = count+1;
    subplot(total_subplots,1,count)
    f_ImageMatrix(squeeze(mean(Time_Freq_Vehicle_subj_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(pink)
    %colorbar
    title('TF Vehicle')

    ValidTimeCue = v_time(~isnan(Time_Freq_Cue_clust.(clusters{cluster})(1,1,:)));
    ValidTimeVehicle = v_time(~isnan(Time_Freq_Vehicle_clust.(clusters{cluster})(1,1,:)));
    
    ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
    start_sample = ceil((ValidTime(1)-v_time(1))/s_tstep)+1;
    end_sample = floor((ValidTime(end)-v_time(1))/s_tstep);
    
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_nonParametricTest_Subj('Beta',...
        ' ',' ','Odor D','Vehicle',...
        v_Beta_Cue_subj.(clusters{cluster}),...
        v_Beta_Vehicle_subj.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample);
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_nonParametricTest_Subj('Spindle',...
        ' ',' ','Odor D','Vehicle',...
        v_Spindle_Cue_subj.(clusters{cluster}),...
        v_Spindle_Vehicle_subj.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample);
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_nonParametricTest_Subj('Theta',...
        ' ',' ','Odor D','Vehicle',...
        v_Theta_Cue_subj.(clusters{cluster}),...
        v_Theta_Vehicle_subj.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample);
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_nonParametricTest_Subj('Delta',...
        'Time(sec)',' ','Odor D','Vehicle',...
        v_Delta_Cue_subj.(clusters{cluster}),...
        v_Delta_Vehicle_subj.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample);
    
    
    set(gcf,'position',[81,35,700,926])
    
    l = suptitle(strcat(filesOdor(subj).name(1:6),'_',clusters{cluster}));
    set(l, 'Interpreter', 'none')
    
    saveas(gcf,strcat(filesOdor(subj).name(1:6),'_',clusters{cluster},'.png'))
    close all

end