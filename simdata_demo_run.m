%====================================
% This script is the main file to run STI, aSTI, Fc-STI, Fc-aSTI on
% simulation data and evaluate the reconstruction results of these models.
%=========================================

%% load data
root_path='./data/simulation';
folder='SNR30_Orinum11_Maxrot30';
mkdir(fullfile(root_path,folder,'output'));   %Create directory for output files
outpath=fullfile(root_path,folder,'output');
% load ground truth
MMS_data=load_nii(fullfile(root_path,'Ground_truth','sim_MMS.nii.gz'));
MMS_GT=MMS_data.img;
MSA_data=load_nii(fullfile(root_path,'Ground_truth','sim_MSA.nii.gz'));
MSA_GT=MSA_data.img;
PEV_data=load_nii(fullfile(root_path,'Ground_truth','sim_pev.nii.gz'));
PEV_GT=PEV_data.img;
% load data for compute STI
phase_data=load_nii(fullfile(root_path,folder,'phi.nii.gz'));   %tissue phase data
phase_tissue=phase_data.img;
H_data=load(fullfile(root_path,folder,'H.mat'));         % Direction of B0 in the subject frame of reference
H_Matrix=H_data.H;
mask_data=load_nii(fullfile(root_path,folder,'mask.nii.gz'));   %brain mask
mask=double(mask_data.img);
WM_mask_data=load_nii(fullfile(root_path,folder,'WM_mask.nii.gz'));   %White matter mask
WM_mask=double(WM_mask_data.img);
FiberDir_data=load_nii(fullfile(root_path,folder,'FiberDir.nii.gz'));   %fiber directions relative to B0: sin^2(theta), used for Fc-aSTI*
FiberDir=FiberDir_data.img;

%% STI model
[chi_tensor_STI,MMS_STI,MSA_STI,PEV_STI]=STI_inverse(phase_tissue,H_Matrix,mask,outpath);
%% aSTI model
[chi_tensor_aSTI,MMS_aSTI,MSA_aSTI,PEV_aSTI]=aSTI_inverse(phase_tissue,H_Matrix,mask,outpath);
%% Fc-STI model
[MMS_FcSTI1,MSA_FcSTI1,PEV_FcSTI1,Fres_FcSTI1]= Fc_STI_iteration(phase_tissue,H_Matrix,mask,WM_mask,1,outpath);
%% Fc-aSTI model
[MMS_FcaSTI1,MSA_FcaSTI1,PEV_FcaSTI1,Fres_FcaSTI1]= Fc_aSTI_iteration(phase_tissue,H_Matrix,mask,WM_mask,1,outpath);
%% assessment
fprintf('MAE of MMS reconstructed by STI: %0.4f\n',compute_mae(MMS_STI,MMS_GT,mask));
fprintf('MAE of MMS reconstructed by aSTI: %0.4f\n',compute_mae(MMS_aSTI,MMS_GT,mask));
fprintf('MAE of MMS reconstructed by Fc-STI*: %0.4f\n',compute_mae(MMS_FcSTI1,MMS_GT,mask));
fprintf('MAE of MMS reconstructed by Fc-aSTI*: %0.4f\n',compute_mae(MMS_FcaSTI1,MMS_GT,mask));

disp('=====================');
fprintf('MAE of MSA reconstructed by STI: %0.4f\n',compute_mae(MSA_STI,MSA_GT,WM_mask));
fprintf('MAE of MSA reconstructed by aSTI: %0.4f\n',compute_mae(MSA_aSTI,MSA_GT,WM_mask));
fprintf('MAE of MSA reconstructed by Fc-STI*: %0.4f\n',compute_mae(MSA_FcSTI1,MSA_GT,WM_mask));
fprintf('MAE of MSA reconstructed by Fc-aSTI*: %0.4f\n',compute_mae(MSA_FcaSTI1,MSA_GT,WM_mask));

[~,AE_mean_STI,~]=compute_AE(PEV_GT,PEV_STI,WM_mask);
[~,AE_mean_aSTI,~]=compute_AE(PEV_GT,PEV_aSTI,WM_mask);
[~,AE_mean_FcSTI1,~]=compute_AE(PEV_GT,PEV_FcSTI1,WM_mask);
[~,AE_mean_FcaSTI1,~]=compute_AE(PEV_GT,PEV_FcaSTI1,WM_mask);

disp('===================');
fprintf('AE of PEV reconstructed by STI: %4.2f\n', AE_mean_STI);
fprintf('AE of PEV reconstructed by aSTI: %4.2f\n', AE_mean_aSTI);
fprintf('AE of PEV reconstructed by Fc-STI*: %4.2f\n', AE_mean_FcSTI1);
fprintf('AE of PEV reconstructed by Fc-aSTI*: %4.2f\n', AE_mean_FcaSTI1);
