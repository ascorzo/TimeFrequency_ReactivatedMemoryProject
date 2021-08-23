addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Odor Night
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

% %files = dir(strcat(filepath,'*.set'));
% filesOdor = dir(strcat(filepath,'*_Odor.mat'));
% filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
% p_clustersOfInterest
% clusters = fieldnames(Clust);


% DNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_DNight');

% for subj = 1:numel(filesOdor)
    
%     v_cycles =  DNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
%     [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);
    
%     s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
%     s_Lastepoch = s_MaxNumStimsIdx;
    
%     load(strcat(filepath,filesOdor(subj).name));
%     load(strcat(filepath,filesVehicle(subj).name));
    
%     cfg                         = [];
%     cfg.trials                  = s_Firstepoch:s_Lastepoch;
%     Time_Freq_Odor_Selected     = ft_selectdata(cfg, Time_Freq_Odor);
%     Time_Freq_Vehi_Selected     = ft_selectdata(cfg, Time_Freq_Vehicle);
    
%     for freq = 1:numel(Time_Freq_Odor_Selected.freq)
%         for chan = 1:numel(Time_Freq_Odor_Selected.label)

%             disp(strcat('Sujeto: ',num2str(subj),{' '},...
%                 'Frequencia: ',num2str(freq),{' '},...
%                 'Channel: ',num2str(chan)))

%             for t = 1:numel(Time_Freq_Odor_Selected.time)

%                 Xt = 1:size(Time_Freq_Odor_Selected.powspctrm,1);

                
%                 OdorEvol    = squeeze(Time_Freq_Odor_Selected.powspctrm(:,chan,freq,t));
%                 VehicleEvol = squeeze(Time_Freq_Vehi_Selected.powspctrm(:,chan,freq,t));
                
%                 mdl_Odor = fitlm(Xt,OdorEvol, 'linear');
%                 mdl_Vehicle  = fitlm(Xt,VehicleEvol, 'linear');
                
%                 % Get slope of regression
%                 slope_DNight_Odor(chan,freq,t) = (mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
%                     (numel(Xt) - 1);
                
%                 slope_DNight_Vehicle(chan,freq,t) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
%                     (numel(Xt) - 1);
%             end
%         end
%     end
    
%     cfg                         = [];
%     cfg.trials                  ='all';
%     cfg.avgoverrpt              ='yes';
%     Time_Freq_Odor_Slope{subj}  = ft_selectdata(cfg, Time_Freq_Odor);
%     Time_Freq_Vehi_Slope{subj}  = ft_selectdata(cfg, Time_Freq_Vehicle);

%     % Replace
%     Time_Freq_Odor_Slope{subj}.powspctrm  = slope_DNight_Odor;
%     Time_Freq_Vehi_Slope{subj}.powspctrm  = slope_DNight_Vehicle;

%     clear slope_DNight_Odor slope_DNight_Vehicle

%     savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

%     save(strcat(savepath,'TF_DNight_TFSlope'),...
%         'Time_Freq_Odor_Slope','Time_Freq_Vehi_Slope')
    
% end



% clear Time_Freq_Odor_Slope Time_Freq_Vehi_Slope 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

MNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_MNight');


for subj = 1:numel(filesOdor)
        
    v_cycles =  MNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
    [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);
    
    s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
    s_Lastepoch = s_MaxNumStimsIdx;
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
    
    cfg                         = [];
    cfg.trials                  = s_Firstepoch:s_Lastepoch;
    Time_Freq_Odor_Selected     = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_Selected     = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    for freq = 1:numel(Time_Freq_Odor_Selected.freq)
        for chan = 1:numel(Time_Freq_Odor_Selected.label)
            disp(strcat('Sujeto: ',num2str(subj),{' '},...
                'Frequencia: ',num2str(freq),{' '},...
                'Channel: ',num2str(chan)))
            for t = 1:numel(Time_Freq_Odor_Selected.time)
                Xt = 1:size(Time_Freq_Odor_Selected.powspctrm,1);

                
                OdorEvol    = squeeze(Time_Freq_Odor_Selected.powspctrm(:,chan,freq,t));
                VehicleEvol = squeeze(Time_Freq_Vehi_Selected.powspctrm(:,chan,freq,t));
                
                mdl_Odor = fitlm(Xt,OdorEvol, 'linear');
                mdl_Vehicle  = fitlm(Xt,VehicleEvol, 'linear');
                
                % Get slope of regression
                slope_DNight_Odor(chan,freq,t) = (mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
                    (numel(Xt) - 1);
                
                slope_DNight_Vehicle(chan,freq,t) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
                    (numel(Xt) - 1);
            end
        end
    end
    
    cfg                         = [];
    cfg.trials                  ='all';
    cfg.avgoverrpt              ='yes';
    Time_Freq_Odor_Slope{subj}  = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_Slope{subj}  = ft_selectdata(cfg, Time_Freq_Vehicle);

    % Replace
    Time_Freq_Odor_Slope{subj}.powspctrm  = slope_DNight_Odor;
    Time_Freq_Vehi_Slope{subj}.powspctrm  = slope_DNight_Vehicle;

    clear slope_DNight_Odor slope_DNight_Vehicle

    savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

    save(strcat(savepath,'TF_MNight_TFSlope'),...
        'Time_Freq_Odor_Slope','Time_Freq_Vehi_Slope')
    
end

