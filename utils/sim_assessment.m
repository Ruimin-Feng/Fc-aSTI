function []=sim_assessment(GTpath,folder)
%=========================
% This function calculates quantitative performance metrics of the results 
% from the four models on simulation data
% Inputs:
% GTpath: directory of the ground truth maps
% folder: directory of the reconstructed maps

%% MMS assessment
MMS_GT_data=load_nii(fullfile(GTpath,'sim_MMS.nii'));
MMS_GT=MMS_GT_data.img;

MMS_symm_data=load_nii(fullfile(folder,'MMS_STI.nii'));
MMS_symm=MMS_symm_data.img;

MMS_symm_Fs_data=load_nii(fullfile(folder,'MMS_Fc-STI.nii'));
MMS_symm_Fs=MMS_symm_Fs_data.img;

MMS_asymm_data=load_nii(fullfile(folder,'MMS_aSTI.nii'));
MMS_asymm=MMS_asymm_data.img;

MMS_asymm_Fs_data=load_nii(fullfile(folder,'MMS_Fc-aSTI.nii'));
MMS_asymm_Fs=MMS_asymm_Fs_data.img;

fprintf('RMSE of MMS reconstructed by STI: %4.2f\n',compute_rmse(MMS_symm,MMS_GT));
fprintf('RMSE of MMS reconstructed by Fc-STI: %4.2f\n',compute_rmse(MMS_symm_Fs,MMS_GT));
fprintf('RMSE of MMS reconstructed by aSTI: %4.2f\n',compute_rmse(MMS_asymm,MMS_GT));
fprintf('RMSE of MMS reconstructed by Fc-aSTI: %4.2f\n',compute_rmse(MMS_asymm_Fs,MMS_GT));

fprintf('SSIM of MMS reconstructed by STI: %4.2f\n',compute_ssim(MMS_symm,MMS_GT));
fprintf('SSIM of MMS reconstructed by Fc-STI: %4.2f\n',compute_ssim(MMS_symm_Fs,MMS_GT));
fprintf('SSIM of MMS reconstructed by aSTI: %4.2f\n',compute_ssim(MMS_asymm,MMS_GT));
fprintf('SSIM of MMS reconstructed by Fc-aSTI: %4.2f\n',compute_ssim(MMS_asymm_Fs,MMS_GT));


%% MSA assessment

MSA_GT_data=load_nii(fullfile(GTpath,'sim_MSA.nii'));
MSA_GT=MSA_GT_data.img;

MSA_symm_data=load_nii(fullfile(folder,'MSA_STI.nii'));
MSA_symm=MSA_symm_data.img;
MSA_symm_Fs_data=load_nii(fullfile(folder,'MSA_Fc-STI.nii'));
MSA_symm_Fs=MSA_symm_Fs_data.img;
MSA_asymm_data=load_nii(fullfile(folder,'MSA_aSTI.nii'));
MSA_asymm=MSA_asymm_data.img;
MSA_asymm_Fs_data=load_nii(fullfile(folder,'MSA_Fc-aSTI.nii'));
MSA_asymm_Fs=MSA_asymm_Fs_data.img;

disp('=====================');
fprintf('RMSE of MSA reconstructed by STI: %4.2f\n',compute_rmse(MSA_symm,MSA_GT));
fprintf('RMSE of MSA reconstructed by Fc-STI: %4.2f\n',compute_rmse(MSA_symm_Fs,MSA_GT));
fprintf('RMSE of MSA reconstructed by aSTI: %4.2f\n',compute_rmse(MSA_asymm,MSA_GT));
fprintf('RMSE of MSA reconstructed by Fc-aSTI: %4.2f\n',compute_rmse(MSA_asymm_Fs,MSA_GT));

fprintf('SSIM of MSA reconstructed by STI: %4.2f\n',compute_ssim(MSA_symm,MSA_GT));
fprintf('SSIM of MSA reconstructed by Fc-STI: %4.2f\n',compute_ssim(MSA_symm_Fs,MSA_GT));
fprintf('SSIM of MSA reconstructed by aSTI: %4.2f\n',compute_ssim(MSA_asymm,MSA_GT));
fprintf('SSIM of MSA reconstructed by Fc-aSTI: %4.2f\n',compute_ssim(MSA_asymm_Fs,MSA_GT));

%% PEV assesssment
PEV_GT_data=load_nii(fullfile(GTpath,'sim_PEV.nii'));
PEV_GT=PEV_GT_data.img;
mask_data=load_nii(fullfile(folder,'mask.nii'));
mask=mask_data.img;

PEV_symm_data=load_nii(fullfile(folder,'PEV_STI.nii'));
PEV_symm=PEV_symm_data.img;

PEV_symm_Fs_data=load_nii(fullfile(folder,'PEV_Fc-STI.nii'));
PEV_symm_Fs=PEV_symm_Fs_data.img;

PEV_asymm_data=load_nii(fullfile(folder,'PEV_aSTI.nii'));
PEV_asymm=PEV_asymm_data.img;

PEV_asymm_Fs_data=load_nii(fullfile(folder,'PEV_Fc-aSTI.nii'));
PEV_asymm_Fs=PEV_asymm_Fs_data.img;

[~,AE_mean_symm,~]=compute_AE(PEV_GT,PEV_symm,mask);
[~,AE_mean_symm_Fs,~]=compute_AE(PEV_GT,PEV_symm_Fs,mask);
[~,AE_mean_asymm,~]=compute_AE(PEV_GT,PEV_asymm,mask);
[~,AE_mean_asymm_Fs,~]=compute_AE(PEV_GT,PEV_asymm_Fs,mask);

disp('===================');
fprintf('AE of PEV reconstructed by STI: %4.2f\n', AE_mean_symm);
fprintf('AE of PEV reconstructed by Fc-STI: %4.2f\n', AE_mean_symm_Fs);
fprintf('AE of PEV reconstructed by aSTI: %4.2f\n', AE_mean_asymm);
fprintf('AE of PEV reconstructed by Fc-aSTI: %4.2f\n', AE_mean_asymm_Fs);

%% offset assessment
Fc_GT_data=load(fullfile(GTpath,'sim_Fs.mat'));
Fc_GT=Fc_GT_data.sim_Fs;

Fc_symm_Fs_data=load_nii(fullfile(folder,'Fc-STI.nii'));
Fc_symm_Fs=Fc_symm_Fs_data.img;
Fc_asymm_Fs_data=load_nii(fullfile(folder,'Fc-aSTI.nii'));
Fc_asymm_Fs=Fc_asymm_Fs_data.img;

disp('=================');
fprintf('RMSE of offset reconstructed by Fc-STI: %4.2f\n',compute_rmse(Fc_symm_Fs(:,:,:,7),Fc_GT));
fprintf('RMSE of offset reconstructed by Fc-aSTI: %4.2f\n',compute_rmse(Fc_asymm_Fs(:,:,:,10),Fc_GT));

fprintf('SSIM of MSA reconstructed by Fc-STI: %4.2f\n',compute_ssim(Fc_symm_Fs(:,:,:,7),Fc_GT));
fprintf('SSIM of MSA reconstructed by Fc-aSTI: %4.2f\n',compute_ssim(Fc_asymm_Fs(:,:,:,10),Fc_GT));




end



