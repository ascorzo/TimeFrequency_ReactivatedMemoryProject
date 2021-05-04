

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1/'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

s_tstep = 0.005; % try with 0.005
s_fstep = 0.005; % 0.005

p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:23%numel(filesOdor)

    %----------------------------------------------------------------------
    %           Create Diff Matrix 
    %----------------------------------------------------------------------
    for cluster = 1:numel(clusters)
        
        % For Declarative Associated Odor Night
        Time_Freq_Diff_D_clust.(clusters{cluster})(subj,:,:) = ...
            Time_Freq_OdorD_clust.(clusters{cluster})(subj,:,:)-...
            Time_Freq_VehicleD_clust.(clusters{cluster})(subj,:,:);
        
        % For Motor Associated Odor Night
        Time_Freq_Diff_M_clust.(clusters{cluster})(subj,:,:) = ...
            Time_Freq_OdorM_clust.(clusters{cluster})(subj,:,:)-...
            Time_Freq_VehicleM_clust.(clusters{cluster})(subj,:,:); 
    end
end



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculation of power change in specific bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Freq = Time_Freq_Odor;
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

SpindleBand     = [10 18];
DeltaBand       = [0.5 4];
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


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
for cluster = 1:numel(clusters)
   
    % For Declarative Associated Odor Night
    v_Spindle_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    
    % For Motor Associated Odor Night
    
    v_Spindle_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,BetaIdx,:),2));
    
end



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cluster = 1:numel(clusters)
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = Time_Freq_Odor.time;
    total_subplots_row      = 2;
    total_subplots_column   = 1;
    count                   = 1;
    frequencies             = Time_Freq_Odor.freq;
    v_xlim                  = [-0.5,3];
    
    TF_Diff_D = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster}),1));
    TF_Diff_M = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster}),1));
    
    bottom  = min([min(TF_Diff_D(:)),min(TF_Diff_M(:))]);
    top     = max([max(TF_Diff_D(:)),max(TF_Diff_M(:))]);
    
    start_sample = 1;
    end_sample = length(time_parcial);
    

    %----------------------------------------------------------------------
    % Plot TF Difference
    %----------------------------------------------------------------------
    
    % -- D Night ----
    tf(count) = subplot(total_subplots_row,total_subplots_column,count,...
        'PlotBoxAspectRatio',[1,0.396226399912024,0.396226399912024]);
%         'Position',[0.13,0.69,0.366055226824458,0.209574269707028],...
         
    f_ImageMatrix(TF_Diff_D,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)','Position',[-6.977238662891239,-3.349991562537163,1.000000000000014])
    colormap(parula)
    caxis manual
    caxis([bottom top]*0.7);
    title('D Night')
    
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         'Position',[0.53,0.69,0.359635108481264,0.209574269707028]);
    f_ImageMatrix(TF_Diff_M,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom top]*0.7);
    title('M Night')
    
    colorbar('Position',...
        [0.92915844819818,0.399715504978663,0.013642340756455,0.500711237553343]);
    
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    
%     count = count+1;
%     tf(count) = subplot(total_subplots_row,total_subplots_column,count);%,...
%         %'Position',[0.13,0.12,0.368027613412229,0.215735294117647]);
%     f_nonParametricTest('',...
%         'Time(sec)',' ','D Night','M Night',...
%         v_Delta_Diff_D.(clusters{cluster}),...
%         v_Delta_Diff_M.(clusters{cluster}),...
%         time_parcial,'-r','-b',start_sample,end_sample);
% %     legend('Position',...
% %         [0.900542877170138,0.256921324754893,0.093688361625934,0.046230439780277]);
    xlim(v_xlim)
    ylim([5 20])
    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x','y')
%     set(gcf,'position',[190,275,1014,703])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    %saveas(gcf,strcat(clusters{cluster},'.png'))

end