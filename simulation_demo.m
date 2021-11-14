%====================================
% This script is a demo of simulation experiments in the paper
% An improved asymmetric susceptibility tensor imaging model with frequency offset correction
% (Feng et al)
%--------------------------------------------
%It includes:
% (1) generate the ground truth susceptibility tensor elements, mean magnetic susceptibility (MMS), 
% magnetic susceptibility anisotropy (MSA), principal eigenvector (PEV), and the simulation phase data with 
% 12 different head orientations with a certain level of Gaussian noise.
% (2) run four different reconstruction models on the simulation data
% (3) calculate quantitative evaluation metrics (RMSE, SSIM, and AE) of the
% results obtained by the four models 
% 
%-------------------------------------------
% Note that NIfTI toolbox should be installed and added to the path to
% ensure .nii files can be successfully read and saved.
%=========================================
addpath('utils');
%% First step
fpath='data/simulation';
load(fullfile(fpath,'dti_tensor.mat'));
load(fullfile(fpath,'dti_mask.mat'));
load(fullfile(fpath,'sim_Fs.mat'));
%  adjusted values of DTI elements to the range of susceptibility tensor values.
init_stitsr=zeros(size(dtitsr));
for i=1:6
    temp=dtitsr(:,:,:,i);
    init_stitsr(:,:,:,i)=((temp-min(temp(:)))/(max(temp(:))-min(temp(:)))*0.5-0.25).*mask;
end
% impose cylindrical symmetry constraint
[Chi_vec1,Chi_eig_init,~,~]=Cal_eig(init_stitsr,mask,'symm');
Chi_eig1=(Chi_eig_init(:,:,:,1)+Chi_eig_init(:,:,:,2))/2;
VVT1(:,:,:,1)=Chi_vec1(:,:,:,1,1).^2;
VVT1(:,:,:,2)=Chi_vec1(:,:,:,1,1).*Chi_vec1(:,:,:,2,1);
VVT1(:,:,:,3)=Chi_vec1(:,:,:,1,1).*Chi_vec1(:,:,:,3,1);
VVT1(:,:,:,4)=Chi_vec1(:,:,:,2,1).*Chi_vec1(:,:,:,2,1);
VVT1(:,:,:,5)=Chi_vec1(:,:,:,2,1).*Chi_vec1(:,:,:,3,1);
VVT1(:,:,:,6)=Chi_vec1(:,:,:,3,1).*Chi_vec1(:,:,:,3,1);
VVT2(:,:,:,1)=Chi_vec1(:,:,:,1,2).^2;
VVT2(:,:,:,2)=Chi_vec1(:,:,:,1,2).*Chi_vec1(:,:,:,2,2);
VVT2(:,:,:,3)=Chi_vec1(:,:,:,1,2).*Chi_vec1(:,:,:,3,2);
VVT2(:,:,:,4)=Chi_vec1(:,:,:,2,2).*Chi_vec1(:,:,:,2,2);
VVT2(:,:,:,5)=Chi_vec1(:,:,:,2,2).*Chi_vec1(:,:,:,3,2);
VVT2(:,:,:,6)=Chi_vec1(:,:,:,3,2).*Chi_vec1(:,:,:,3,2);
VVT3(:,:,:,1)=Chi_vec1(:,:,:,1,3).^2;
VVT3(:,:,:,2)=Chi_vec1(:,:,:,1,3).*Chi_vec1(:,:,:,2,3);
VVT3(:,:,:,3)=Chi_vec1(:,:,:,1,3).*Chi_vec1(:,:,:,3,3);
VVT3(:,:,:,4)=Chi_vec1(:,:,:,2,3).*Chi_vec1(:,:,:,2,3);
VVT3(:,:,:,5)=Chi_vec1(:,:,:,2,3).*Chi_vec1(:,:,:,3,3);
VVT3(:,:,:,6)=Chi_vec1(:,:,:,3,3).*Chi_vec1(:,:,:,3,3);
sim_stitsr=(Chi_eig1.*VVT1+Chi_eig1.*VVT2+Chi_eig_init(:,:,:,3).*VVT3).*mask;
save_nii(make_nii(sim_stitsr.*mask),fullfile(fpath,'sim_stitsr.nii'));
save_nii(make_nii(mask),fullfile(fpath,'mask.nii'));
[Chi_vec,Chi_eig,MMS,MSA]=Cal_eig(sim_stitsr,mask,'symm');
eigenVec3=Chi_vec(:,:,:,:,3);
cMSA=abs(eigenVec3).*repmat(MSA,[1 1 1 3]);
save_nii(make_nii(cMSA),fullfile(fpath,'sim_cMSA.nii'));
save_nii(make_nii(MMS),fullfile(fpath,'sim_MMS.nii'));
save_nii(make_nii(MSA),fullfile(fpath,'sim_MSA.nii'));
save_nii(make_nii(Chi_eig),fullfile(fpath,'sim_eig.nii'));
save_nii(make_nii(eigenVec3),fullfile(fpath,'sim_PEV.nii'));
save_nii(make_nii(abs(eigenVec3)),fullfile(fpath,'sim_absPEV.nii'));

%simulate frequency data
Ndir=12;  %12 head orientations
NoiseStd=0.05;  % determine a Gaussian noise level
outpath=fullfile(fpath,strcat(num2str(NoiseStd*100),'%_Noise'));
mkdir(outpath);
save_nii(make_nii(mask),fullfile(outpath,'mask.nii'));
H=zeros(Ndir,3);
for ori=1:Ndir
    theta_x=rand(1)*2*pi-pi;
    theta_y=rand(1)*2*pi-pi;
    Mrot=[1,0,0;0,cos(theta_x),-sin(theta_x);0,sin(theta_x),cos(theta_x)]*[cos(theta_y),0,sin(theta_y);0,1,0;-sin(theta_y),0,cos(theta_y)];
    temp=Mrot*[0,0,1]';
    H(ori,:)=temp';
end
perturb=STI_forward(sim_stitsr,H);
phi=zeros(size(perturb));
for i=1:size(phi,4)
    temp=perturb(:,:,:,i);  
    Noise=NoiseStd*max(temp(:))*randn(size(temp));
    phi(:,:,:,i)=(perturb(:,:,:,i)+sim_Fs+Noise).*mask;
end
save_nii(make_nii(phi),fullfile(outpath,'phi.nii'));
save(fullfile(outpath,'phi.mat'),'phi');
save(fullfile(outpath,'H.mat'),'H');

%% second step
run_models(phi, H, mask, outpath);

%% third step
sim_assessment(fpath,outpath);
