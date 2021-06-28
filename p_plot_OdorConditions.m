%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1'))
addpath('./Scripts_Wilc')

tf = [];
p_clustersOfInterest
clusters = fieldnames(Clust);

%clusters = {'All'};

for cluster = 1:numel(clusters)
    
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = Time_Freq_Odor.time;
    total_subplots_row      = 6;
    total_subplots_column   = 1;
    count                   = 1;
    frequencies             = Time_Freq_Odor.freq;
    v_xlim                  = [-5 25];
    
    TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
    TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD(:)),...
        min(TF_OdorM(:))])*0.6;
    top_TF     = max([max(TF_OdorD(:)),...
        max(TF_OdorM(:))])*0.6;
    
    bottom_Delta  = min([min(v_Delta_OdorD.(clusters{cluster})(:)),...
        min(v_Delta_OdorM.(clusters{cluster})(:))])*0.2;
    
    top_Delta     = max([max(v_Delta_OdorD.(clusters{cluster})(:)),...
        max(v_Delta_OdorM.(clusters{cluster})(:))])*0.2;
    
    start_sample = 1;
    end_sample = length(time_parcial);
    

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
    
    % -- D Night ----
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)','position',[-6.77975582249371,-6.043966945036731,1])
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor D')
    xlabel('')
    
    
    % -- M Night ----  
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor M')
    xlabel('')
    a = colorbar('Position',...
    [0.923139158749658,0.659827213822894,0.011812297561022,0.23110151187905]);
    a.Label.String = 'Zscore';
    a.Label.Position = [0.365367938365255,0.14086177304198,0];
    a.Label.Rotation = 0;
      
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),'-',num2str(SpindleBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Odor M',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot Theta Power
    %----------------------------------------------------------------------

    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),'-',num2str(ThetaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Odor M',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    ylabel('Zscore','position',[-6.95692119341007,-0.171399629483115,-0.999999999999986])

    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------

    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),'-',num2str(DeltaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Odor M',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot SW Power
    %----------------------------------------------------------------------
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),'-',num2str(SWBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Odor M',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
    xlim(v_xlim)
    legend('Position',...
        [0.751703720334021,0.037836151274367,0.158396948020877,0.040631750283952]);
    
    xticks([-5 0 5 10 15 20 25])
    

    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    set(gcf,'position',[533.8,41.8,524,740.8])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    savepath = 'G:\Mi unidad\2021\AnalysisTemp\Figures\OdorDvsM_TF\';
    saveas(gcf,strcat(savepath,clusters{cluster},'.png'))

end