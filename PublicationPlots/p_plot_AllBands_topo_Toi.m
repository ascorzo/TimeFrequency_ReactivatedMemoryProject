% plot topo of toi all bands

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);


%=========================================================================
% In order to be more organized with the results, and not have multiple 
% TF figures for each cluster of channels, I would like to check the
% topographical map of the power for the defined frequency bands in the
% critical points where I can see some significance
%--------------------------------------------------------------------------

% The moments of interest are:
% SO: [-1 1],[9 12],[9.5 10.5],[11 13],[15 20]
% Delta: [-0.5 0.5],[3.5 4.5],[7 8],[9 11],[15 20]
% Theta: [-0.5 0.5],[-4.5 -3.5],[7 8],[10 12]
% Spindles: [2 3]

Toi.SW          = [-1 1;9 12;9.5 10.5;11 13;15 20];
Toi.Delta       = [-0.5 0.5;3.5 4.5;7 8;9 11;15 20];
Toi.Theta       = [-0.5 0.5;-4.5 -3.5;7 8;10 12];
Toi.Spindle     = [2 3];

bands = fieldnames(Toi);

SpindleBand     = [12 16];
DeltaBand       = [1 4];
ThetaBand       = [4 8];
SWBand          = [0.5 2];


for subj = 1:numel(filesOdor)
    
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
    for band = bands
        
        % select frequency band
        
        
        
    end
end
