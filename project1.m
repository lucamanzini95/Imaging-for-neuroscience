clear 
close all
clc

% Define base paths for the 'old' and 'young' datasets
basePathOld = 'dataset/old';
basePathYoung = 'dataset/young';

% Load the Schaefer segmentation and white matter mask (assuming these are common across subjects)
parcelAtlas = load_untouch_nii('dataset/Schaefer_segmentation.nii.gz');
wmMask = load_untouch_nii('dataset/MNI_WM_mask.nii.gz');
wmMaskLogical = logical(wmMask.img);


%% TASK 1

% Define the number of parcels
numParcels = max(parcelAtlas.img(:));

% Process both groups
processGroup(basePathOld, 'old', numParcels, parcelAtlas, wmMaskLogical);
processGroup(basePathYoung, 'young', numParcels, parcelAtlas, wmMaskLogical);

% Function definition for processing each group
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

%% TASK 2

load dataset/labels.mat
old_avg=load('old_averageSignals3D.mat');
old_avg=old_avg.groupAverageSignals;
old_avg=permute(old_avg,[2 1 3]); % change of dimension to compute the correlation correclty

young_avg=load('young_averageSignals3D.mat');
young_avg=young_avg.groupAverageSignals;
young_avg=permute(young_avg,[2 1 3]);

[n_volumes, ~ , n_subj]=size(old_avg);

%% OLD (TASK 2)

for subj=1:n_subj
    FC_old_static(:,:,subj)=corr(old_avg(:,:,subj)); % correlation (Pearson)
end
   
zFC_old_static=atanh(FC_old_static); % Fisher's z-transform

figure
for subj=1:n_subj
    subplot(2,5,subj)
    imagesc(zFC_old_static(:,:,subj))
    title(['Subject n: ',num2str(subj)])
    %set(gca, 'XTick', 1:n_parcel, 'XTickLabel', labels, 'XTickLabelRotation', 90, 'FontSize', 3)
    %set(gca, 'YTick', 1:n_parcel, 'YTickLabel', labels, 'FontSize', 3)
end
colormap jet
sgtitle('Static FC old')

%% YOUNG (TASK 2)


for subj=1:n_subj
    FC_young_static(:,:,subj)=corr(young_avg(:,:,subj));
end

zFC_young_static=atanh(FC_young_static);

figure
for subj=1:n_subj
    subplot(2,5,subj)
    imagesc(zFC_young_static(:,:,subj))
    title(['Subject n: ',num2str(subj)])
    %set(gca, 'XTick', 1:n_parcel, 'XTickLabel', labels, 'XTickLabelRotation', 90, 'FontSize', 3)
    %set(gca, 'YTick', 1:n_parcel, 'YTickLabel', labels, 'FontSize', 3)
end
colormap jet
sgtitle('Static FC young')

%% TASKS 3 & 4

%% OLD (TASKS 3 & 4)

i=1;
window=30;
for subj=1:n_subj   
    for time=1:2:n_volumes-window
        windowed_signal=old_avg(time:time+window,:,subj); % to get the windowed signal
        FC_old_dynamic(:,:,subj,i)=corr(windowed_signal);
        i=i+1;
    end
    i=1;
end

zFC_old_dynamic=atanh(FC_old_dynamic);
numParcels=size(zFC_old_dynamic,1);

%% YOUNG (TASKS 3 & 4)

i=1;
window=30;
for subj=1:n_subj   
    for time=1:2:n_volumes-window
        windowed_signal=young_avg(time:time+window,:,subj); % to get the windozwed signal
        FC_young_dynamic(:,:,subj,i)=corr(windowed_signal);
        i=i+1;
    end
    i=1;
end                                                                                                                                             

zFC_young_dynamic=atanh(FC_young_dynamic);

%% TASK 5: Characterization of the old population

% Step 5a: Pack the upper triangular values of the zFC matrices into a matrix for clustering
numWindows = size(zFC_old_dynamic, 4);
upperTriIndices = find(triu(ones(numParcels), 1));  % Indices of the upper triangular part
old_packed_FC = zeros(numWindows * n_subj, length(upperTriIndices));  % Initialize packed matrix

% Fill the packed matrix with the upper triangular values
index = 1;
for subj = 1:n_subj
    for window = 1:numWindows
        tempMatrix = zFC_old_dynamic(:, :, subj, window);
        old_packed_FC(index, :) = tempMatrix(upperTriIndices);
        index = index + 1;
    end
end

% Step 5b: Perform clustering using hierarchical clustering and determine the optimal number of clusters
Z = linkage(old_packed_FC, 'ward');
optimalClusters = 2:10;  % Test a range of clusters
silhouetteValues = zeros(length(optimalClusters), 1);

for k = optimalClusters
    clusterIdx = cluster(Z, 'maxclust', k);
    silhouetteValues(k-1) = mean(silhouette(old_packed_FC, clusterIdx));
end

% Determine the optimal number of clusters
[~, bestK] = max(silhouetteValues);
K_Old = optimalClusters(bestK);

% Perform clustering with the optimal number of clusters
clusterIdx = cluster(Z, 'maxclust', K_Old);

% Step 5c: Compute fractional occupancy
fractionalOccupancy = zeros(n_subj, K_Old);

for k = 1:K_Old
    for subj = 1:n_subj
        % Calculate the proportion of time windows assigned to cluster k for each subject
        fractionalOccupancy(subj, k) = sum(clusterIdx(((subj-1)*numWindows+1):(subj*numWindows)) == k) / numWindows;
    end
end

% Step 5d: Characterize each state
A_Old = zeros(numParcels, numParcels, K_Old);

for k = 1:K_Old
    % Calculate the centroid of each cluster
    centroid = mean(old_packed_FC(clusterIdx == k, :), 1);
    % Reconstruct the FC matrix from the upper triangular values
    centroidMatrix = zeros(numParcels, numParcels);
    centroidMatrix(upperTriIndices) = centroid;
    centroidMatrix = centroidMatrix + centroidMatrix';  % Symmetrize the matrix

    % Step 5d(i): Reduce the density
    threshold = prctile(centroidMatrix(upperTriIndices), 75);  % 75th percentile threshold
    sparseMatrix = centroidMatrix >= threshold;
    A_Old(:, :, k) = sparseMatrix;

    % Step 5d(ii): Compute the Laplacian matrix
    D = diag(sum(A_Old(:, :, k), 2));  % Degree matrix
    L = D - A_Old(:, :, k);  % Laplacian matrix

    % Step 5d(iii): Compute the principal eigenvector of the Laplacian matrix
    [eigVecs, eigVals] = eig(L);
    eigVals = diag(eigVals);  % Extract eigenvalues
    [~, minIdx] = min(eigVals);  % Index of the smallest eigenvalue
    principalEigVec = eigVecs(:, minIdx);  % Principal eigenvector

    % Save the gradient vector
    Gradient_Old(:, k) = principalEigVec;
end

% Save the results
save('A_Old.mat', 'A_Old');
save('Gradient_Old.mat', 'Gradient_Old');
save('fractionalOccupancy_Old.mat', 'fractionalOccupancy');

% Visualize the results
figure;
for k = 1:K_Old
    subplot(2, ceil(K_Old/2), k);
    imagesc(A_Old(:, :, k));
    title(['State ' num2str(k)]);
    axis square;
    colorbar;
end
sgtitle('Sparse Matrices for Old Population States');

%%
%% TASK 6: Characterization of the young population

% Step 5a: Pack the upper triangular values of the zFC matrices into a matrix for clustering
numWindows = size(zFC_young_dynamic, 4);
upperTriIndices = find(triu(ones(numParcels), 1));  % Indices of the upper triangular part
young_packed_FC = zeros(numWindows * n_subj, length(upperTriIndices));  % Initialize packed matrix

% Fill the packed matrix with the upper triangular values
index = 1;
for subj = 1:n_subj
    for window = 1:numWindows
        tempMatrix = zFC_young_dynamic(:, :, subj, window);
        young_packed_FC(index, :) = tempMatrix(upperTriIndices);
        index = index + 1;
    end
end

% Step 5b: Perform clustering using hierarchical clustering and determine the optimal number of clusters
Z = linkage(young_packed_FC, 'ward');
optimalClusters = 2:10;  % Test a range of clusters
silhouetteValues = zeros(length(optimalClusters), 1);

for k = optimalClusters
    clusterIdx = cluster(Z, 'maxclust', k);
    silhouetteValues(k-1) = mean(silhouette(young_packed_FC, clusterIdx));
end

% Determine the optimal number of clusters
[~, bestK_young] = max(silhouetteValues);
K_Young = optimalClusters(bestK_young);

% Perform clustering with the optimal number of clusters
clusterIdx = cluster(Z, 'maxclust', K_Young);

% Step 5c: Compute fractional occupancy
fractionalOccupancy = zeros(n_subj, K_Young);

for k = 1:K_Young
    for subj = 1:n_subj
        % Calculate the proportion of time windows assigned to cluster k for each subject
        fractionalOccupancy(subj, k) = sum(clusterIdx(((subj-1)*numWindows+1):(subj*numWindows)) == k) / numWindows;
    end
end

% Step 5d: Characterize each state
A_Young = zeros(numParcels, numParcels, K_Young);

for k = 1:K_Young
    % Calculate the centroid of each cluster
    centroid = mean(young_packed_FC(clusterIdx == k, :), 1);
    % Reconstruct the FC matrix from the upper triangular values
    centroidMatrix = zeros(numParcels, numParcels);
    centroidMatrix(upperTriIndices) = centroid;
    centroidMatrix = centroidMatrix + centroidMatrix';  % Symmetrize the matrix

    % Step 5d(i): Reduce the density
    threshold = prctile(centroidMatrix(upperTriIndices), 75);  % 75th percentile threshold
    sparseMatrix = centroidMatrix >= threshold;
    A_Young(:, :, k) = sparseMatrix;

    % Step 5d(ii): Compute the Laplacian matrix
    D = diag(sum(A_Young(:, :, k), 2));  % Degree matrix
    L = D - A_Young(:, :, k);  % Laplacian matrix

    % Step 5d(iii): Compute the principal eigenvector of the Laplacian matrix
    [eigVecs, eigVals] = eig(L);
    eigVals = diag(eigVals);  % Extract eigenvalues
    [~, minIdx] = min(eigVals);  % Index of the smallest eigenvalue
    principalEigVec = eigVecs(:, minIdx);  % Principal eigenvector

    % Save the gradient vector
    Gradient_Young(:, k) = principalEigVec;
end

% Save the results
save('A_Young.mat', 'A_Young');
save('Gradient_Young.mat', 'Gradient_Young');
save('fractionalOccupancy_Young.mat', 'fractionalOccupancy');

% Visualize the results
figure;
for k = 1:K_Young
    subplot(2, ceil(K_Young/2), k);
    imagesc(A_Young(:, :, k));
    title(['State ' num2str(k)]);
    axis square;
    colorbar;
end
sgtitle('Sparse Matrices for Young Population States');

%% Task 7
% Load the data for young population
load('A_Young.mat');
load('Gradient_Young.mat');
load('fractionalOccupancy_Young.mat');

% Assuming you have K_Young defined similarly as K_Old
K_Young = size(A_Young, 3);  % Number of states for the young population

% Compare K_Old and K_Young
if K_Old == K_Young
    disp('The number of clusters (states) is the same for both old and young populations.');
else
    disp('The number of clusters (states) differs between old and young populations.');
end

% Combine gradients from both populations
allGradients = [Gradient_Old, Gradient_Young];

% Compute Spearman's correlation between all gradients
Dist_All = corr(allGradients, 'Type', 'Spearman');

% Visualize the correlation matrix
figure;
imagesc(Dist_All);
colorbar;
title('Spearman Correlation Matrix of Gradients (Old and Young)');
xlabel('States (Old + Young)');
ylabel('States (Old + Young)');

% Extract the correlation between old and young states
correlationOldYoung = Dist_All(1:K_Old, K_Old+1:end);

% Find the minimum correlation value for each old state
[minCorrelation, minIndex] = min(correlationOldYoung, [], 2);

% Identify the old state with the least similarity to any young state
[~, leastSimilarStateIndex] = min(minCorrelation);
leastSimilarState = leastSimilarStateIndex;

% Display the results
disp(['The state of the old population that least resembles the youth states is state ' num2str(leastSimilarState)]);
disp(['The minimum correlation value for this state is ' num2str(minCorrelation(leastSimilarState))]);

% Provide a metric to quantify the decision
metric = minCorrelation(leastSimilarState);
disp(['Metric (Minimum Correlation) for the least similar state: ' num2str(metric)]);
