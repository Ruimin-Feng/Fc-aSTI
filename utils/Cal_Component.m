function []=Cal_Component(path)
%=======================
% This function calculates eigenvalues, eigenvectors, MMS, and MSA from the
% susceptibility tensor elements reconstructed by the four models at the same time.
% Inputs:
% path: directory where the susceptibility tensors reconstructed by the
% four models located
%=======================
path0=pwd;
cd(path);
symm_data=load_nii('STI.nii');
symm=symm_data.img;
symm_Fs_data=load_nii('Fc-STI.nii');
symm_Fs=symm_Fs_data.img;
asymm_data=load_nii('aSTI_symmpart.nii');
asymm=asymm_data.img;
asymm_Fs_data=load_nii('Fc-aSTI_symmpart.nii');
asymm_Fs=asymm_Fs_data.img;
mask_data=load_nii(fullfile('mask.nii'));
mask=mask_data.img;

[Vec_symm,Eig_symm,MMS_symm,MSA_symm]=Cal_eig(symm,mask,'symm');
[Vec_symm_Fs,Eig_symm_Fs,MMS_symm_Fs,MSA_symm_Fs]=Cal_eig(symm_Fs,mask,'symm');
[Vec_asymm,Eig_asymm,MMS_asymm,MSA_asymm]=Cal_eig(asymm,mask,'asymm');
[Vec_asymm_Fs,Eig_asymm_Fs,MMS_asymm_Fs,MSA_asymm_Fs]=Cal_eig(asymm_Fs,mask,'asymm');

save_nii(make_nii(MMS_symm),'MMS_STI.nii');
save_nii(make_nii(MSA_symm),'MSA_STI.nii');
save_nii(make_nii(Eig_symm),'Eig_STI.nii');
save_nii(make_nii(Vec_symm),'Vec_STI.nii');

save_nii(make_nii(MMS_symm_Fs),'MMS_Fc-STI.nii');
save_nii(make_nii(MSA_symm_Fs),'MSA_Fc-STI.nii');
save_nii(make_nii(Eig_symm_Fs),'Eig_Fc-STI.nii');
save_nii(make_nii(Vec_symm_Fs),'Vec_Fc-STI.nii');

save_nii(make_nii(MMS_asymm),'MMS_aSTI.nii');
save_nii(make_nii(MSA_asymm),'MSA_aSTI.nii');
save_nii(make_nii(Eig_asymm),'Eig_aSTI.nii');
save_nii(make_nii(Vec_asymm),'Vec_aSTI.nii');

save_nii(make_nii(MMS_asymm_Fs),'MMS_Fc-aSTI.nii');
save_nii(make_nii(MSA_asymm_Fs),'MSA_Fc-aSTI.nii');
save_nii(make_nii(Eig_asymm_Fs),'Eig_Fc-aSTI.nii');
save_nii(make_nii(Vec_asymm_Fs),'Vec_Fc-aSTI.nii');
cd(path0);
end



