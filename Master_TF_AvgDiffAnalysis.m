addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off


filepath_D = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';
filesOdor_D = dir(strcat(filepath_D,'*_Odor.mat'));
filesVehicle_D = dir(strcat(filepath_D,'*_Vehicle.mat'));


filepath_M = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';
filesOdor_M = dir(strcat(filepath_M,'*_Odor.mat'));
filesVehicle_M = dir(strcat(filepath_M,'*_Vehicle.mat'));

p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor_D)
    
    disp(strcat('Sujeto: ',num2str(subj)))

    %----------------------------------------------------------------------
    % For Declarative Associated Odor Night
    %----------------------------------------------------------------------

    %---------- Load data---------------------
    load(strcat(filepath_D,filesOdor_D(subj).name));
    load(strcat(filepath_D,filesVehicle_D(subj).name));

    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';

    Time_Freq_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_baseline2);

    Time_Freq_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_baseline2);
    
    
    %----------- Analysis by clusters --------------------------

    for cluster = 1:numel(clusters)

        %------ Combine Mean of channels by clusters -----------------------
        chans = intersect(Clust.(clusters{cluster}),Time_Freq_Odor.label);
        cfg = [];
        cfg.channel = chans;
        cfg.avgoverchan = 'yes';
        Time_Freq_Odor_clust = ...
            ft_selectdata(cfg,Time_Freq_Odor);

        Time_Freq_Vehicle_clust = ...
            ft_selectdata(cfg,Time_Freq_Vehicle);

        %---------- Difference between Odor and Vehicle ---------------------
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.operation    = '(x1-x2)';
        Time_Freq_Diff_D.(clusters{cluster}){subj} = ...
            ft_math(cfg, Time_Freq_Odor_clust, Time_Freq_Vehicle_clust);

        Time_Freq_Diff_D_clust.(clusters{cluster})(subj,:,:) = ...
            Time_Freq_Diff_D.(clusters{cluster}){subj}.powspctrm;
    end


    %----------------------------------------------------------------------
    % For Motor Associated Odor Night
    %----------------------------------------------------------------------

    %---------- Load data---------------------
    load(strcat(filepath_M,filesOdor_M(subj).name));
    load(strcat(filepath_M,filesVehicle_M(subj).name));

    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';

    Time_Freq_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_baseline2);

    Time_Freq_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_baseline2);
    
    
    %----------- Analysis by clusters --------------------------

    for cluster = 1:numel(clusters)

        %------ Combine Mean of channels by clusters -----------------------
        chans = intersect(Clust.(clusters{cluster}),Time_Freq_Odor.label);
        cfg = [];
        cfg.channel = chans;
        cfg.avgoverchan = 'yes';
        Time_Freq_Odor_clust = ...
            ft_selectdata(cfg,Time_Freq_Odor);

        Time_Freq_Vehicle_clust = ...
            ft_selectdata(cfg,Time_Freq_Vehicle);

        %---------- Difference between Odor and Vehicle ---------------------
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.operation    = '(x1-x2)';
        Time_Freq_Diff_M.(clusters{cluster}){subj} = ...
            ft_math(cfg, Time_Freq_Odor_clust, Time_Freq_Vehicle_clust);  

        Time_Freq_Diff_M_clust.(clusters{cluster})(subjh,:,:) = ...
            Time_Freq_Diff_M.(clusters{cluster}){subj}.powspctrm;
    end 
end

%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       Calculation of power change in specific bands
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Freq = Time_Freq_Odor;


SpindleBand     = [10 18];
DeltaBand       = [0.5 2];
ThetaBand       = [4 8];
BetaBand        = [18 30];
time            = Time_Freq.time ;
subjects        = 1:23;%numel(filesOdor);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
BetaIdx = find(Time_Freq.freq>=BetaBand(1) &...
    Time_Freq.freq<=BetaBand(2));


for cluster = 1:numel(clusters)
    % -----------------For D Night ----------------------------------
    v_Spindle_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Diff_D.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

    % -----------------For M Night----------------------------------
    v_Spindle_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Diff_M.(clusters{cluster}) = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];

for cluster = 1:numel(clusters)
    
    figure
    y_lims          = [];
    x_lims_parcial  = [-0.5 3];
    time_parcial    = Time_Freq.time;
    total_subplots  = 4;
    count           = 1;
    frequencies     = Time_Freq.freq;
    
    TF_D_Night = squeeze(mean(Time_Freq_Diff_D_clust.(clusters{cluster}),1));
    TF_M_Night = squeeze(mean(Time_Freq_Diff_M_clust.(clusters{cluster}),1));   
    
    bottom  = min([min(TF_D_Night(:)),min(TF_M_Night(:))])*0.7;
    top     = max([max(TF_D_Night(:)),max(TF_M_Night(:))])*0.7;
    
    %----------------------------------------------------------------------
    % For Declarative Associated Odor Night
    %----------------------------------------------------------------------

    tf(count) = subplot(total_subplots,1,count);
    f_ImageMatrix(TF_D_Night,time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(parula)
    caxis manual
    caxis([bottom top]);
    %colorbar
    title('D Night')
    
 
    count = count+1;
    tf(count) = subplot(total_subplots,1,count);
    f_ImageMatrix(TF_M_Night,time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(parula)
    caxis manual
    caxis([bottom top]);
    %colorbar
    title('M Night')
    
    linkaxes(tf,'x')
    
    %set(gcf,'position',[81,35,700,926])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    %saveas(gcf,strcat(clusters{cluster},'Final_OdorD.png'))
    
    
    v_time = Time_Freq.time ;
    ValidTime_DNight = v_time(~isnan(Time_Freq_Diff_D_clust.(clusters{cluster})(1,1,:)));
    ValidTime_MNight = v_time(~isnan(Time_Freq_Diff_M_clust.(clusters{cluster})(1,1,:)));
    
    ValidTime = intersect(ValidTime_DNight,ValidTime_MNight);
%     start_sample = ceil((ValidTime(1)-v_time(1))/s_tstep)+1;
%     end_sample = floor((ValidTime(end)-v_time(1))/s_tstep);

    start_sample = 1;
    end_sample = length(v_time);

%     count = count+1;
%     subplot(total_subplots,1,count)
%     f_nonParametricTest('Beta',...
%         ' ',' ','D Night','M Night',...
%         v_Beta_Diff_D.(clusters{cluster}),...
%         v_Beta_Diff_M.(clusters{cluster}),...
%         v_time,'-r','-b',start_sample,end_sample);
%     
    count = count+1;
    tf(count) = subplot(total_subplots,1,count);
    f_nonParametricTest('Spindle',...
        ' ',' ','D Night','M Night',...
        v_Spindle_Diff_D.(clusters{cluster}),...
        v_Spindle_Diff_M.(clusters{cluster}),...
        v_time,'-r','-b',start_sample,end_sample);
    xlim(x_lims_parcial)
    
    count = count+1;
    tf(count) = subplot(total_subplots,1,count);
    f_nonParametricTest('Theta',...
        ' ',' ','D Night','M Night',...
        v_Theta_Diff_D.(clusters{cluster}),...
        v_Theta_Diff_M.(clusters{cluster}),...
        v_time,'-r','-b',start_sample,end_sample);
    xlim(x_lims_parcial)
    
%     count = count+1;
%     tf(count) = subplot(total_subplots,1,count);
%     f_nonParametricTest('Delta',...
%         'Time(sec)',' ','D Night','M Night',...
%         v_Delta_Diff_D.(clusters{cluster}),...
%         v_Delta_Diff_M.(clusters{cluster}),...
%         v_time,'-r','-b',start_sample,end_sample);
%     xlim(x_lims_parcial)
   
    
    linkaxes(tf,'x')
end
