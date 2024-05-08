# Imaging-for-neuroscience

The goal of this homework is to investigate the relationship between brain states of two different groups of individuals as defined by dynamic functional connectivity obtained from resting state functional MRI data.

### Dataset

DATASET
Dataset included:
old = folder containing data of 10 old healthy controls  
young = folder containing data of 10 young healthy controls  


For both old and young controls, each subfolder contains:  
preproc_rs_fMRI.nii.gz = preprocessed fMRI data (2x2x2 mm3, TR=1 s, 400 volumes) mapped into the MNI152 symmetric atlas  
brain_mask.nii.gz = individual fMRI brain mask mapped into the MNI152 symmetric atlas


Data also included:
Schaefer_segmentation.nii.gz = Schaefer segmentation (100 parcels - 7 RSNs) of the MNI152 symmetric atlas
Schaefer2018_100Parcels_7Networks_order.txt = file containing the labels of the parcels. For example, parcel 1 (7Networks_LH_Vis_1) is referring to a
parcel of the Visual Network (RSN) of the Left Hemisphere labels.mat = vector with the labels of the parcels
MNI_WM_mask.nii.gz = white matter mask of the MNI152 symmetric atlas