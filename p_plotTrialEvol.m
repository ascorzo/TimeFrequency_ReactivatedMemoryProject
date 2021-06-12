

load('G:\Mi unidad\2021\AnalysisTemp\firstSec_TFcount_DNight.mat');
load('G:\Mi unidad\2021\AnalysisTemp\firstSec_TFcount_MNight.mat');

p_clustersOfInterest
clusters = fieldnames(Clust);


for cluster = 1:numel(clusters)
    Time_Freq_OdorD_evol.(clusters{cluster}) = Time_Freq_OdorD_evol.(clusters{cluster})(:,1:7,:);
    Time_Freq_OdorM_evol.(clusters{cluster}) = Time_Freq_OdorM_evol.(clusters{cluster})(:,1:7,:);
    Time_Freq_VehicleD_evol.(clusters{cluster}) = Time_Freq_VehicleD_evol.(clusters{cluster})(:,1:7,:);
    Time_Freq_VehicleM_evol.(clusters{cluster}) = Time_Freq_VehicleM_evol.(clusters{cluster})(:,1:7,:);
end

Xt = 1:size(Time_Freq_OdorD_evol.central,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% D Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cluster = 1:numel(clusters) 
    figure
    
    
    nplots = 5;
    count = 1;
    
    % SW plot
    data1 = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(1,:,:))';
    data2 = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(1,:,:))';
    
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('SW')
    
    % Delta plot
    data1 = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(2,:,:))';
    data2 = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(2,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Delta')
    
    % Theta plot
    data1 = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(3,:,:))';
    data2 = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(3,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Theta')
    
    % Slow Spindles plot
    data1 = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(4,:,:))';
    data2 = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(4,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Slow Spindles')
    
    % Fast Spindles plot
    data1 = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(5,:,:))';
    data2 = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(5,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Fast Spindles')
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    set(gcf,'position',[488,41.800000000000004,560,740.8000000000001])
end
close all
%%

%for Theta

for cluster = 1:numel(clusters)
    
    figure
    
    OdorEvol = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(3,:,:));
    VehicleEvol = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(3,:,:));
    
    subplot(3,1,1)
    bar(mean(OdorEvol,2))
    hold on
    errorbar(1:size(OdorEvol,1),mean(OdorEvol,2),std(OdorEvol,[],2)/size(OdorEvol,2))
    title('Odor D')
        
    subplot(3,1,2)
    bar(mean(VehicleEvol,2))
    hold on
    errorbar(1:size(VehicleEvol,1),mean(VehicleEvol,2),std(VehicleEvol,[],2)/size(VehicleEvol,2))
    title('Vehicle')
    
    subplot(3,1,3)
    Difference = OdorEvol-...
        VehicleEvol;
    
    mean_Diff   = mean(Difference,2);
    errbar      = std(Difference,[],2);
    sembar = errbar/size(Difference,2);
    
    bar(mean_Diff)
    hold on
    errorbar(1:length(errbar),mean_Diff,sembar)
    title('Odor D  -  Vehicle')
    xlabel('Trials')
    l = sgtitle(strcat(clusters{cluster},{' '},'Theta z-score'));
    set(l, 'Interpreter', 'none')
    
    set(gcf,'position',[488,41.800000000000004,560,740.8000000000001])
    
    saveas(gcf,strcat((clusters{cluster}),'_DNight','.png'))
end
close all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xt = 1:size(Time_Freq_OdorM_evol.central,2);

for cluster = 1:numel(clusters) 
    figure
    
    
    nplots = 5;
    count = 1;
    
    % SW plot
    data1 = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(1,:,:))';
    data2 = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(1,:,:))';
    
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('SW')
    
    % Delta plot
    data1 = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(2,:,:))';
    data2 = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(2,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Delta')
    
    % Theta plot
    data1 = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(3,:,:))';
    data2 = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(3,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Theta')
    
    % Slow Spindles plot
    data1 = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(4,:,:))';
    data2 = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(4,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Slow Spindles')
    
    % Fast Spindles plot
    data1 = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(5,:,:))';
    data2 = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(5,:,:))';
    
    count = count+1;
    subplot(nplots,1,count)
    p1 = shadedErrorBar(Xt,mean(data1,1),...
        std(data1,[],1)/sqrt(size(data1,1)),'lineprops', '-r');
    hold on
    p2 = shadedErrorBar(Xt,mean(data2,1),...
        std(data2,[],1)/sqrt(size(data1,1)),'lineprops', '-k');
    title('Fast Spindles')
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    set(gcf,'position',[488,41.800000000000004,560,740.8000000000001])
    
    
end
close all
%%

%for Theta

for cluster = 1:numel(clusters)
    
    figure
    
    OdorEvol = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(3,:,:));
    VehicleEvol = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(3,:,:));
    
    subplot(3,1,1)
    bar(mean(OdorEvol,2))
    hold on
    errorbar(1:size(OdorEvol,1),mean(OdorEvol,2),std(OdorEvol,[],2)/size(OdorEvol,2))
    title('Odor M')
        
    subplot(3,1,2)
    bar(mean(VehicleEvol,2))
    hold on
    errorbar(1:size(VehicleEvol,1),mean(VehicleEvol,2),std(VehicleEvol,[],2)/size(VehicleEvol,2))
    title('Vehicle')
    
    subplot(3,1,3)
    Difference = OdorEvol-...
        VehicleEvol;
    
    mean_Diff   = mean(Difference,2);
    errbar      = std(Difference,[],2);
    sembar = errbar/size(Difference,2);
    
    bar(mean_Diff)
    hold on
    errorbar(1:length(errbar),mean_Diff,sembar)
    title('Odor M  -  Vehicle')
    xlabel('Trials')
    
    l = sgtitle(strcat(clusters{cluster},{' '},'Theta z-score'));
    set(l, 'Interpreter', 'none')
    
    set(gcf,'position',[488,41.800000000000004,560,740.8000000000001])
    
    saveas(gcf,strcat((clusters{cluster}),'_MNight','.png'))
end
close all

%% Implementar esto por sujeto para ver la regresión del cambio de theta

for cluster = 1:numel(clusters)
    
    for subj = 1:23
        
        OdorEvol    = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(3,:,subj));
        VehicleEvol = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(3,:,subj));
        DiffEvol    = OdorEvol-VehicleEvol;
        
        X = 1:size(Time_Freq_OdorD_evol.(clusters{cluster}),2);
        mdl_Odor = fitlm(X,OdorEvol, 'linear');
        mdl_Vehicle  = fitlm(X,VehicleEvol, 'linear');
        mdl_Diff  = fitlm(X,DiffEvol, 'linear');
        
        %
        %     figure
        %     plot(OdorEvol,'r')
        %     hold on
        plot(mdl_Odor.Fitted,'k')
        %
        %
        %     hold on
        %     plot(VehicleEvol,'k')
        %     hold on
        %     plot(mdl_Vehicle.Fitted,'k')
        
        %     plot(mdl_Diff.Fitted,'b')
        hold on
        
        
        % Get slope of regression
        slope_DNight_Odor(subj) = ( mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_DNight_Vehicle(subj) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_DNight_Diff(subj) = ( mdl_Diff.Fitted(end) - mdl_Diff.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        
    end
    
    OdorEvol = squeeze(Time_Freq_OdorD_evol.(clusters{cluster})(3,:,:));
    VehicleEvol = squeeze(Time_Freq_VehicleD_evol.(clusters{cluster})(3,:,:));
    
    DiffEvol = mean(OdorEvol-VehicleEvol,2);
    mdl_Diff  = fitlm(X,mean(OdorEvol,2), 'linear');
    plot(mdl_Diff.Fitted,'r')
    close all
    
    
    %% Implementar esto por sujeto para ver la regresión del cambio de theta
    
    
    
    
    figure
    for subj = 1:23
        
        OdorEvol    = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(3,:,subj));
        VehicleEvol = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(3,:,subj));
        DiffEvol    = OdorEvol-VehicleEvol;
        
        X = Xt;
        
        mdl_Odor = fitlm(X,OdorEvol, 'linear');
        mdl_Vehicle  = fitlm(X,VehicleEvol, 'linear');
        mdl_Diff  = fitlm(X,DiffEvol, 'linear');
        
        plot(mdl_Odor.Fitted,'k')
        hold on
        
        
        % Get slope of regression
        slope_MNight_Odor(subj) = ( mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Vehicle(subj) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Diff(subj) = ( mdl_Diff.Fitted(end) - mdl_Diff.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
    end
    
    OdorEvol = squeeze(Time_Freq_OdorM_evol.(clusters{cluster})(3,:,:));
    VehicleEvol = squeeze(Time_Freq_VehicleM_evol.(clusters{cluster})(3,:,:));
    
    DiffEvol = mean(OdorEvol-VehicleEvol,2);
    mdl_Diff  = fitlm(X,mean(OdorEvol,2), 'linear');
    plot(mdl_Diff.Fitted,'r')
    
    %%
    
    figure
    
    pl(1) = subplot(1,3,1);
    Odor_slopes = [slope_DNight_Odor' slope_MNight_Odor'];
    boxplot(Odor_slopes)
    [h,p_OdorSlope] = ttest(slope_DNight_Odor,slope_MNight_Odor);
    title({'Odor',strcat('p =',num2str(p_OdorSlope))})
    xticklabels({'D Night','M Night'})
    
    pl(2) = subplot(1,3,2);
    Vehicle_slopes = [slope_DNight_Vehicle' slope_MNight_Vehicle'];
    boxplot(Vehicle_slopes)
    [h,p_VehicleSlope] = ttest(slope_DNight_Vehicle,slope_MNight_Vehicle);
    title({'Vehicle',strcat('p =',num2str(p_VehicleSlope))})
    xticklabels({'D Night','M Night'})
    
    pl(3) = subplot(1,3,3);
    Diff_slopes = [slope_DNight_Diff' slope_MNight_Diff'];
    boxplot(Diff_slopes)
    [h,p_DiffSlope] = ttest(slope_DNight_Diff,slope_MNight_Diff);
    title({'Odor - Vehicle',strcat('p =',num2str(p_DiffSlope))})
    xticklabels({'D Night','M Night'})
    linkaxes(pl,'y')
    
    saveas(gcf,strcat('Theta_evolution_slope_',(clusters{cluster}),'.png'))
    
    
    % Slope odor vs Vehicle
    
    figure
    subplot(1,2,1)
    OvsV_slopes = [slope_DNight_Odor' slope_DNight_Vehicle'];
    boxplot(OvsV_slopes)
    [h,p_OvsVSlope] = ttest(slope_DNight_Odor,slope_DNight_Vehicle);
    xticklabels({'Odor D','Vehicle'})
    title({'D Night',strcat('p =',num2str(p_OvsVSlope))})
    
    
    subplot(1,2,2)
    OvsV_slopes = [slope_MNight_Odor' slope_MNight_Vehicle'];
    boxplot(OvsV_slopes)
    [h,p_OvsVSlope] = ttest(slope_MNight_Odor,slope_MNight_Vehicle);
    xticklabels({'Odor M','Vehicle'})
    title({'M Night',strcat('p =',num2str(p_OvsVSlope))})
    
    saveas(gcf,strcat('Theta_evolution_slope_OdorvsVehicle',(clusters{cluster}),'.png'))
    
end
close all


