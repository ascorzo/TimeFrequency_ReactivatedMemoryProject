addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

s_tstep = 0.1; % try with 0.005
s_fstep = 0.15; % 0.005

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Organize files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/DNight/';


%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));

for subj = 1:numel(filesOdor)
    
  
    % disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   
    if ~isempty(find(isnan(Time_Freq_Odor.powspctrm)))
        disp(filesOdor(subj).name)
    end
    
    if ~isempty(find(isnan(Time_Freq_Vehicle.powspctrm)))
        disp(filesVehicle(subj).name)
    end   

    disp('Done')
     
end


%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/MNight/';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));

for subj = 1:numel(filesOdor)

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   
       
    if ~isempty(find(isnan(Time_Freq_Odor.powspctrm)))
        disp(filesOdor(subj).name)
    end
    
    if ~isempty(find(isnan(Time_Freq_Vehicle.powspctrm)))
        disp(filesVehicle(subj).name)
    end 
     
    disp('Done')
end
