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