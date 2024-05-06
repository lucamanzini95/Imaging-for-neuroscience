clear all
close all
clc

%setting paths
Utils_path= '/Users/luca/Documents/MATLAB/NeuroImag/Utils'; %insert the path where you stored the Utils/NIFTI_toolbox folder
addpath(genpath(Utils_path));%adding to matlab path

% Define base paths for the 'old' and 'young' datasets
basePathOld = 'dataset/old';
basePathYoung = 'dataset/young';

% Load the Schaefer segmentation and white matter mask 
parcelAtlas = load_untouch_nii('dataset/Schaefer_segmentation.nii.gz');
wmMask = load_untouch_nii('dataset/MNI_WM_mask.nii.gz');
wmMaskLogical = logical(wmMask.img);
%% Task 1: for each healthy control (HC) extract the average fMRI signals for each parcel of the atlas
% exclude voxel outside the individual mask
%exclude voxel within the white matter task

% Define the number of parcels
numParcels = max(parcelAtlas.img(:));

% Processing function for each group, now including parcelAtlas and wmMaskLogical as arguments
function processGroup(basePath, groupName, numParcels, parcelAtlas, wmMaskLogical)
    % List all subject folders in the base path
    subjects = dir(basePath);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));  % Filter to keep only directories
    
    % Initialize matrix to store all subjects' average signals
    groupAverageSignals = zeros(numParcels, length(subjects));
    
    % Loop over each subject
    for i = 1:length(subjects)
        % Load the subject's fMRI data and brain mask
        fmriDataPath = fullfile(basePath, subjects(i).name, 'preproc_rs_fMRI.nii.gz');
        brainMaskPath = fullfile(basePath, subjects(i).name, 'brain_mask.nii.gz');
        
        fmriData = load_untouch_nii(fmriDataPath);
        brainMask = load_untouch_nii(brainMaskPath);
        
        % Loop over each parcel to compute the average fMRI signal
        for parcelID = 1:numParcels
            % Create a mask for the current parcel
            parcelMask = parcelAtlas.img == parcelID;
            
            % Exclude voxels outside the brain mask and within the white matter mask
            finalMask = parcelMask & logical(brainMask.img) & ~wmMaskLogical;
            
            % Extract the fMRI data for the current parcel using the final mask
            parcelData = fmriData.img(finalMask);
            
            % Compute the average signal only if there are any voxels in the mask
            if sum(finalMask(:)) > 0
                groupAverageSignals(parcelID, i) = mean(parcelData);
            else
                groupAverageSignals(parcelID, i) = NaN; % Assign NaN if the parcel has no valid voxels
            end
        end
    end
    
    % Save the group's average signals to a file
    save([groupName '_averageSignals.mat'], 'groupAverageSignals');
end

% Process both groups
processGroup(basePathOld, 'old', numParcels, parcelAtlas, wmMaskLogical);
processGroup(basePathYoung, 'young', numParcels, parcelAtlas, wmMaskLogical);