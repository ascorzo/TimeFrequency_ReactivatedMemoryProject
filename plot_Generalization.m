addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light

F = Time_Freq_Odor.freq;

for subj = 1:numel(cf_Generalization)
    
    figure
    %F = Time_Freq_Odor.freq;
    mv_plot_2D(cf_Generalization{subj}, 'x', F, 'y', F)
    xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
    title('Frequency generalization using channels-x-times as features')
    colorbar('off')
    colorbar
    caxis([min(cf_Generalization{subj}(:)),...
        max(cf_Generalization{subj}(:))])
    
    saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\Odor_DvsM_Generalization\AroundCero\',...
        'Subj',num2str(subj),'.png'))
    
    close all
end
% 
result_average_Generalization = mv_combine_results(resultGeneralization, 'average');

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(result_average_Generalization.perf{1}, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization using channels x time as features')
colorbar('off')
colorbar
caxis([min(result_average_Generalization.perf{1}(:)),...
        max(result_average_Generalization.perf{1}(:))])

saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\Odor_DvsM_Generalization\AroundCero\',...
        'Avegare.png'))
close all


%% statistics

defaultSubj = load('G:\Mi unidad\2021\AnalysisTemp\RC_512_sleepTF_DN_Odor.mat');

cfg             = [];
cfg.trials      = 1;
cfg.channel     = 'all';
cfg.avgoverchan = 'yes';
defaultSubj     = ft_selectdata(cfg, defaultSubj.Time_Freq_Odor);

filename     = 'G:\Mi unidad\2021\AnalysisTemp\dummyfile.set';
sensors      = ft_read_sens(filename);

sensors.chanpos = sensors.chanpos(1,:);
sensors.chantype = sensors.chantype(1,:);
sensors.chanunit = sensors.chanunit(1,:);
sensors.elecpos = sensors.elecpos(1,:);
sensors.label = sensors.label(1,:);


%%

Generalization_DNight = load('G:\Mi unidad\2021\AnalysisTemp\OdorDvsVehicle_FreqGeneralization_[-1,5].mat');
Generalization_MNight = load('G:\Mi unidad\2021\AnalysisTemp\OdorMvsVehicle_FreqGeneralization_[-1,5].mat');

F = TF_D_Generalization{1}.freq;
% completar para D Night
for subj = 1:numel(Generalization_DNight.cf_Generalization)
    TF_D_Generalization{subj}             = defaultSubj;
    TF_D_Generalization{subj}.time        = TF_D_Generalization{subj}.freq;
    TF_D_Generalization{subj}.powspctrm   = TF_D_Generalization{subj}.powspctrm(:,:,:,1:numel(F));
    TF_D_Generalization{subj}.powspctrm(1,1,:,:)   = Generalization_DNight.cf_Generalization{subj}; 
end

% completar para M Night
for subj = 1:numel(Generalization_MNight.cf_Generalization)
    TF_M_Generalization{subj}             = defaultSubj;
    TF_M_Generalization{subj}.time        = TF_M_Generalization{subj}.freq;
    TF_M_Generalization{subj}.powspctrm   = TF_M_Generalization{subj}.powspctrm(:,:,:,1:numel(F));
    TF_M_Generalization{subj}.powspctrm(1,1,:,:)   = Generalization_MNight.cf_Generalization{subj}; 
end

cfg                     = [];
cfg.latency             = 'all';
cfg.frequency           = 'all';
cfg.channel             = 'all';
cfg.correctm            = 'cluster';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';     % use actvsblT for activation against baseline
cfg.clusterstatistic    = 'maxsum';            % statistic used to decide cluster significance (sum of t-values within a cluster)

cfg.tail                = 0;
cfg.clustertail         = cfg.tail;
cfg.alpha               = 0.025;
cfg.numrandomization    = 1000;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;

cfg_neighb.method       = 'distance';
cfg_neighb.channel      = 'all';
cfg_neighb.elec         = sensors; 
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb);

% Design the statistical contrast
design                  = [];
design(1,:)             = [1:length(TF_D_Generalization) 1:length(TF_M_Generalization)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(TF_D_Generalization)) ones(1,length(TF_M_Generalization))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

stats1                  = ft_freqstatistics(cfg, TF_D_Generalization{:}, TF_M_Generalization{:});



%% plot supersubject

addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light

F = Time_Freq_Odor.freq;

    % Odor D vs Vehicle
load('G:\Mi unidad\2021\AnalysisTemp\OdorDvsVehicle_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor D vs Vehicle')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])
    
saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\OdorDvsVehicleGeneralization\AroundCero\',...
        'SuperSubject.png'))
close all

    % Odor M vs Vehicle
    
load('G:\Mi unidad\2021\AnalysisTemp\OdorMvsVehicle_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor M vs Vehicle')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])
saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\OdorMvsVehicleGeneralization\AroundCero\',...
        'SuperSubject.png'))
close all

% Odor D vs Odor M
load('G:\Mi unidad\2021\AnalysisTemp\OdorDvsOdorM_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor D vs Odor M')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])
saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\Odor_DvsM_Generalization\AroundCero\',...
        'SuperSubject.png'))
close all


%% plot supersubject Above threshold

addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light

F = Time_Freq_Odor.freq;

    % Odor D vs Vehicle
load('G:\Mi unidad\2021\AnalysisTemp\OdorDvsVehicle_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf>=0.54, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor D vs Vehicle')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])
    


    % Odor M vs Vehicle
    
load('G:\Mi unidad\2021\AnalysisTemp\OdorMvsVehicle_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf>=0.54, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor M vs Vehicle')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])


% Odor D vs Odor M
load('G:\Mi unidad\2021\AnalysisTemp\OdorDvsOdorM_FreqGeneralization_[-1,5]_SuperSubj.mat')

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(resultGeneralization.perf>=0.54, 'x', F, 'y', F)
xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
title('Frequency generalization Odor D vs Odor M')
colorbar('off')
colorbar
caxis([min(resultGeneralization.perf(:)),...
        max(resultGeneralization.perf(:))])

