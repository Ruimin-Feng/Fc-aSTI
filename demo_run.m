%====================================
% This script is the main file to run STI, aSTI, Fc-STI, Fc-aSTI
%=========================================

fpath='./data/mouse_brain';
mkdir(fullfile(fpath,'output'));
outpath=fullfile(fpath,'output');
%% load data
mask_data=load(fullfile(fpath,'mask.mat'));
mask=mask_data.msk;
phase_data=load(fullfile(fpath,'Phi_Prohance1.mat'));
phase_tissue=phase_data.Phi_Prohance1/(2*pi*42.575*7*0.008);
H_data=load(fullfile(fpath,'H_Prohance1.mat'));
H_Matrix=H_data.H_Prohance1;    
offset_mask_data=load_nii(fullfile(fpath,'offset_mask.nii.gz'));
offset_mask=offset_mask_data.img;
%% STI model
[chi_tensor_STI,MMS_STI,MSA_STI,PEV_STI]=STI_inverse(phase_tissue,H_Matrix,mask,outpath);
%% aSTI model
[chi_tensor_aSTI,MMS_aSTI,MSA_aSTI,PEV_aSTI]=aSTI_inverse(phase_tissue,H_Matrix,mask,outpath);
%% Fc-STI* model
[MMS_FcSTI1,MSA_FcSTI1,PEV_FcSTI1,Fres_FcSTI1]= Fc_STI_iteration(phase_tissue,H_Matrix,mask,offset_mask,1,outpath);
%% Fc-aSTI* model
[MMS_FcaSTI1,MSA_FcaSTI1,PEV_FcaSTI1,Fres_FcaSTI1]= Fc_aSTI_iteration(phase_tissue,H_Matrix,mask,offset_mask,1,outpath);