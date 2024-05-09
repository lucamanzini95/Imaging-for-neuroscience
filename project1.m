clear all
close all
clc

% Define base paths for the 'old' and 'young' datasets
basePathOld = 'dataset/old';
basePathYoung = 'dataset/young';

% Load the Schaefer segmentation and white matter mask (assuming these are common across subjects)
parcelAtlas = load_untouch_nii('dataset/Schaefer_segmentation.nii.gz');
wmMask = load_untouch_nii('dataset/MNI_WM_mask.nii.gz');
wmMaskLogical = logical(wmMask.img);

% Define the number of parcels
numParcels = max(parcelAtlas.img(:));

%% TASK 1
% Process both groups, now passing additional necessary variables
processGroup(basePathOld, 'old', numParcels, parcelAtlas, wmMaskLogical);
processGroup(basePathYoung, 'young', numParcels, parcelAtlas, wmMaskLogical);

% Function definition for processing each group, now with additional parameters
function processGroup(basePath, groupName, numParcels, parcelAtlas, wmMaskLogical)

    % List all subject folders in the base path
    subjects = dir(basePath);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));  % Filter to keep only directories
    
    % Initialize matrix to store all subjects' average signals
    % Assuming 400 volumes and 10 subjects as a typical setup
    groupAverageSignals = zeros(numParcels, 400, length(subjects));
    
    % Loop over each subject
    for i = 1:length(subjects)
        % Load the subject's fMRI data and brain mask
        fmriDataPath = fullfile(basePath, subjects(i).name, 'preproc_rs_fMRI.nii.gz');
        brainMaskPath = fullfile(basePath, subjects(i).name, 'brain_mask.nii.gz');
        
        fmriData = load_untouch_nii(fmriDataPath);
        brainMask = load_untouch_nii(brainMaskPath);
        
        % Loop over each parcel to compute the average fMRI signal for each volume
        for parcelID = 1:numParcels
            parcelMask = parcelAtlas.img == parcelID;
            finalMask = parcelMask & logical(brainMask.img) & ~wmMaskLogical;
            
            for vol = 1:size(fmriData.img, 4)  % Assuming fourth dimension is time
                parcelData = fmriData.img(:,:,:,vol);
                maskedData = parcelData(finalMask);
                
                if sum(finalMask(:)) > 0
                    groupAverageSignals(parcelID, vol, i) = mean(maskedData);
                else
                    groupAverageSignals(parcelID, vol, i) = NaN; % Assign NaN if no valid voxels
                end
            end
        end
    end
    
    % Save the group's average signals to a file
    save([groupName '_averageSignals3D.mat'], 'groupAverageSignals');
end

%% Task 2
% Load processed data for the old group
load('old_averageSignals3D.mat');  % This loads the variable `groupAverageSignals`
oldData = groupAverageSignals;  % Assuming the data is stored in this variable

% Load processed data for the young group
load('young_averageSignals3D.mat');  % This loads the variable `groupAverageSignals`
youngData = groupAverageSignals;  % Assuming the data is stored in this variable

parcelLabels = readtable('dataset/Schaefer2018_100Parcels_7Networks_order.txt', 'ReadVariableNames', false);
parcelLabels = parcelLabels.Var1;  % Extracting labels from the table

% Compute and visualize functional connectivity for the old group
computeAndVisualizeFC(oldData, parcelLabels, 'Old');

% Compute and visualize functional connectivity for the young group
computeAndVisualizeFC(youngData, parcelLabels, 'Young');

function computeAndVisualizeFC(groupData, groupLabels, groupName)
    numSubjects = size(groupData, 3);
    numParcels = size(groupData, 1);
    numVolumes = size(groupData, 2);
    
    % Preallocate memory for FC matrices and z-transformed FC matrices
    FC_matrices = zeros(numParcels, numParcels, numSubjects);
    zFC_matrices = zeros(numParcels, numParcels, numSubjects);

    % Compute FC for each subject
    for i = 1:numSubjects
        subjectData = squeeze(groupData(:, :, i));  % Time series data for current subject
        FC = corrcoef(subjectData');  % Compute correlation matrix

        % Fisher z-transform of the correlation coefficients
        zFC = atanh(FC);  
        zFC(isinf(zFC)) = 0;  % Handle infinite values which might result from correlations of 1 or -1

        FC_matrices(:, :, i) = FC;
        zFC_matrices(:, :, i) = zFC;
    end

    % Visualization of zFC matrices for each subject in this group
    figure('Name', ['zFC Matrices for ' groupName], 'NumberTitle', 'off');
    for i = 1:numSubjects
        subplot(2, 5, i);
        imagesc(zFC_matrices(:, :, i), [-0.5, 0.5]);  % Scale for visualization
        colorbar;
        title(['Subject ' num2str(i)]);
        xlabel('Parcel Index');
        ylabel('Parcel Index');
        axis square;
    end

    sgtitle([groupName ' Group - zFC Matrices']);  % Super title for the figure

    % Save the FC matrices and zFC matrices to files
    save([groupName '_FC_matrices.mat'], 'FC_matrices');
    save([groupName '_zFC_matrices.mat'], 'zFC_matrices');
end
