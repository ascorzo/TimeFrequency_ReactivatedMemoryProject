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
    total_subplots_row      = 6;
    total_subplots_column   = 1;
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.13,0.69,0.366055226824458,0.209574269707028],...
%         'PlotBoxAspectRatio',[1,0.396226399912024,0.396226399912024]);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)')%,'Position',[-6.977238662891239,-3.349991562537163,1.000000000000014])
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor D')
    
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%          'Position',[0.53,0.69,0.359635108481264,0.209574269707028]);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor M')
    ylabel('Frequency(Hz)')
      
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
    f_nonParametricTest('Spindle',...
        'Time(sec)',' ','Odor D','Odor M',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
%     legend('Position',...
%         [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot Theta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
    f_nonParametricTest('Theta',...
        'Time(sec)',' ','Odor D','Odor M',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
%     legend('Position',...
%         [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot around 2.5 Hz Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
    f_nonParametricTest('Around 2.5Hz',...
        'Time(sec)',' ','Odor D','Odor M',...
        v_Around_2_5_OdorD.(clusters{cluster}),...
        v_Around_2_5_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
%     legend('Position',...
%         [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
    f_nonParametricTest('Theta',...
        'Time(sec)',' ','Odor D','Odor M',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_OdorM.(clusters{cluster}),...
        time_parcial,'-r','-b',start_sample,end_sample);
%     legend('Position',...
%         [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    xticks([-5 0 5 10 15 20 25])

    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    %set(gcf,'position',[190,275,1014,703])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    %saveas(gcf,strcat(clusters{cluster},'.png'))

end