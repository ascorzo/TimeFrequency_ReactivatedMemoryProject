addpath('C:\Users\lanan\Documents\Github\TimeFrequency_ReactivatedMemoryProject\')

p_clustersOfInterest
clusters = fieldnames(Clust);
clusters{length(clusters)+1}  = 'all';

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Clust_ToPlot_DNight.mat')
load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Clust_ToPlot_MNight.mat')

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_MeanAllchans_ToPlot_DNight.mat')
load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_MeanAllchans_ToPlot_MNight.mat')

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Parms')

Time_Freq_OdorD_clust.all = Time_Freq_OdorD_All;
Time_Freq_VehicleD_clust.all = Time_Freq_VehicleD_All;

Time_Freq_OdorM_clust.all = Time_Freq_OdorM_All;
Time_Freq_VehicleM_clust.all = Time_Freq_VehicleM_All;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculation of power change in specific bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Freq = Time_Freq_Parms;
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

SpindleBand     = [12 16];
DeltaBand       = [1 4];
Around_2_5      = [2 8];
ThetaBand       = [4 8];
SWBand          = [0.5 2];
%BetaBand       = [18 30];
time            = Time_Freq.time;
subjects        = 1:23; %numel(filesOdor);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
Around_2_5Idx = find(Time_Freq.freq>=Around_2_5(1) &...
    Time_Freq.freq<=Around_2_5(2));
SWIdx = find(Time_Freq.freq>=SWBand(1) &...
    Time_Freq.freq<=SWBand(2));

%clusters = {'All'};
%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------
for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
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
    v_SW_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
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
    time_parcial            = Time_Freq_Parms.time;
    total_subplots_row      = 6;
    total_subplots_column   = 3;
    count                   = 1;
    frequencies             = Time_Freq_Parms.freq;
    v_xlim                  = [-5 20];
    v_xticks                = [-5 0 5 10 15 20 25];
    
    
    TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
    TF_VehicleD = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster}),1));
    
    TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
    TF_VehicleM = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD(:)),min(TF_VehicleD(:)),...
        min(TF_OdorM(:)),min(TF_VehicleM(:))])*0.6;
    top_TF     = max([max(TF_OdorD(:)),max(TF_VehicleD(:)),...
        max(TF_OdorM(:)),max(TF_VehicleM(:))])*0.6;
    
    max_TF = max(abs(bottom_TF),abs(top_TF));
    
    % -- ylims Spindles
    bottom_Spindles  = min([min(v_Spindle_OdorD.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        min(v_Spindle_OdorM.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Spindles     = max([max(v_Spindle_OdorD.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        max(v_Spindle_OdorM.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Theta
    bottom_Theta  = min([min(v_Theta_OdorD.(clusters{cluster})(:)),...
        min(v_Theta_VehicleD.(clusters{cluster})(:)),...
        min(v_Theta_OdorM.(clusters{cluster})(:)),...
        min(v_Theta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Theta     = max([max(v_Theta_OdorD.(clusters{cluster})(:)),...
        max(v_Theta_VehicleD.(clusters{cluster})(:)),...
        max(v_Theta_OdorM.(clusters{cluster})(:)),...
        max(v_Theta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Delta
    bottom_Delta  = min([min(v_Delta_OdorD.(clusters{cluster})(:)),...
        min(v_Delta_VehicleD.(clusters{cluster})(:)),...
        min(v_Delta_OdorM.(clusters{cluster})(:)),...
        min(v_Delta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Delta     = max([max(v_Delta_OdorD.(clusters{cluster})(:)),...
        max(v_Delta_VehicleD.(clusters{cluster})(:)),...
        max(v_Delta_OdorM.(clusters{cluster})(:)),...
        max(v_Delta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims SW
    bottom_SW  = min([min(v_SW_OdorD.(clusters{cluster})(:)),...
        min(v_SW_VehicleD.(clusters{cluster})(:)),...
        min(v_SW_OdorM.(clusters{cluster})(:)),...
        min(v_SW_VehicleM.(clusters{cluster})(:))])*0.2;
    top_SW     = max([max(v_SW_OdorD.(clusters{cluster})(:)),...
        max(v_SW_VehicleD.(clusters{cluster})(:)),...
        max(v_SW_OdorM.(clusters{cluster})(:)),...
        max(v_SW_VehicleM.(clusters{cluster})(:))])*0.12;
    
    
    start_sample = 1;
    end_sample = length(time_parcial);
    

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
   
    
    % -- M Night ----  
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,2,...
        'Position',[0.5233798195242,0.855206042570128,0.305988515176375,0.089345535800919]);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Odor M')
    xlabel('')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,5,...
        'Position',[0.525840853158327,0.72179636669222,0.305988515176375,0.089345535800919]);
    f_ImageMatrix(TF_VehicleM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Vehicle')
    xlabel('Time (s)')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
        % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,1,...
        'Position',[0.171222313371618,0.855206042570114,0.306218211648892,0.089345535800919]);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    ylabel('Frequency (Hz)','Position',[-6.884988259404048,-3.446376583590895,1])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Odor D')
    xlabel('')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,4,...
        'Position',[0.171452009844135,0.719076938180593,0.305373256767845,0.0923576200811]);
    
    f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Vehicle')
    xlabel('Time (s)')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
%     % -- Both Odors ----  
%     
%     tf(count) = subplot(total_subplots_row,total_subplots_column,3);
%     f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
%     xlim(v_xlim)
%     colormap(parula)
%     caxis manual
%     caxis([bottom_TF top_TF]);
%     title('Odor D')
%     xlabel('')
%     
%     
%     count = count+1;
%     tf(count) = subplot(total_subplots_row,total_subplots_column,6);
%     f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
%     xlim(v_xlim)
%     colormap(parula)
%     caxis manual
%     caxis([bottom_TF top_TF]);
%     title('Odor M')
%     xlabel('')
    a = colorbar('Position',[0.846735996404601,0.713822894168464,0.011321300281366,0.235421162817687]);
    a.Label.String = 'Power (Z)';
    
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,7);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,8);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Spindle_OdorM.(clusters{cluster}),...
        v_Spindle_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,9);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot Theta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,10);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    ylabel('Power (Z)','position',[-9.497365001570415,-0.116481838204675,-1])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,11);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Theta_OdorM.(clusters{cluster}),...
        v_Theta_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,12);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,13);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,14);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Delta_OdorM.(clusters{cluster}),...
        v_Delta_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
   
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,15);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot SW Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,16);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.265927658827671,0.038916064881278,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xticks([-5 0 5 10 15 20 25]) 
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,17);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_SW_OdorM.(clusters{cluster}),...
        v_SW_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.545970189311835,0.041075892258916,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,18);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend('Position',...
        [0.825707679057529,0.039995978652004,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    

    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    set(gcf,'position',[345.8,41.8,975.2,740.8])
    
%     l = sgtitle(clusters{cluster});
%     set(l, 'Interpreter', 'none')
    
    savepath = 'G:\Mi unidad\2021\AnalysisTemp\Figures\TF_Publication\AllConds_fig\';
    saveas(gcf,strcat(savepath,clusters{cluster},'.png'))
    close all
end







