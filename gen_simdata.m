%====================================
% This script is to generate simulation data in the paper
% An improved asymmetric susceptibility tensor imaging model with frequency offset correction
% (Feng et al)
%--------------------------------------------
%-------------------------------------------
% Note that NIfTI toolbox should be installed and added to the path to
% ensure .nii or .nii.gz files can be successfully read and saved.
%=========================================

%% load data
addpath('utils');
fpath='data/simulation';
FA_data=load_nii(fullfile(fpath,'Phantom','dti_FA.nii.gz'));
FA=FA_data.img;
dtipev_data=load_nii(fullfile(fpath,'Phantom','dti_pev.nii.gz'));
dtipev=dtipev_data.img;
mask_data=load_nii(fullfile(fpath,'Phantom','mask.nii.gz'));
mask=double(mask_data.img);
WM_mask_data=load_nii(fullfile(fpath,'Phantom','WM_mask.nii.gz'));
WM_mask=WM_mask_data.img;
DGM_mask_data=load_nii(fullfile(fpath,'Phantom','DGM_mask.nii.gz'));
DGM_mask=DGM_mask_data.img;

%% simulate chi=chi_iso+chi_aniso
mkdir(fullfile(fpath,'Ground_truth'));
[s1,s2,s3]=size(mask);
chi_iso=zeros([s1,s2,s3]);
chi_iso(DGM_mask==1)=0.15;
chi_iso(DGM_mask==2)=0.05;
chi_iso(DGM_mask==3)=0.05;
chi_iso(DGM_mask==4)=0;
chi_iso(DGM_mask==5)=0.03;
chi_iso(WM_mask==1)=-0.05;
save_nii(make_nii(chi_iso),fullfile(fpath,'Ground_truth','chi_iso.nii.gz'));
WM_mask=double(WM_mask);
FAnorm=FA/0.59.*WM_mask;
chi_aniso=0.01*WM_mask.*FAnorm;
save_nii(make_nii(chi_aniso),fullfile(fpath,'Ground_truth','chi_aniso.nii.gz'));
%simulate susceptibility tensor in the subject frame
dtipev=dtipev.*WM_mask;
VVT1(:,:,:,1)=dtipev(:,:,:,1).^2;
VVT1(:,:,:,2)=dtipev(:,:,:,1).*dtipev(:,:,:,2);
VVT1(:,:,:,3)=dtipev(:,:,:,1).*dtipev(:,:,:,3);
VVT1(:,:,:,4)=dtipev(:,:,:,2).*dtipev(:,:,:,2);
VVT1(:,:,:,5)=dtipev(:,:,:,2).*dtipev(:,:,:,3);
VVT1(:,:,:,6)=dtipev(:,:,:,3).*dtipev(:,:,:,3);
IdentityMtx=repmat(reshape([1,0,0,1,0,1],1,1,1,6),s1,s2,s3,1).*WM_mask;
sim_stitsr=(chi_aniso.*VVT1-chi_aniso/2.*(IdentityMtx-VVT1)).*mask;
sim_stitsr(:,:,:,1)=sim_stitsr(:,:,:,1)+chi_iso;
sim_stitsr(:,:,:,4)=sim_stitsr(:,:,:,4)+chi_iso;
sim_stitsr(:,:,:,6)=sim_stitsr(:,:,:,6)+chi_iso;
save_nii(make_nii(sim_stitsr),fullfile(fpath,'Ground_truth','sim_stitsr.nii.gz'));
save_nii(make_nii(dtipev),fullfile(fpath,'Ground_truth','sim_pev.nii.gz'));
save_nii(make_nii(abs(dtipev)),fullfile(fpath,'Ground_truth','sim_abspev.nii.gz'));
save_nii(make_nii(chi_iso),fullfile(fpath,'Ground_truth','sim_MMS.nii.gz'));
save_nii(make_nii(3*chi_aniso/2),fullfile(fpath,'Ground_truth','sim_MSA.nii.gz'));
save_nii(make_nii(mask),fullfile(fpath,'Ground_truth','mask.nii.gz'));

%% simulate frequency data
Ndir=11;  %number of head orientations (Ndir>=6)
SNR=30;  % Gaussian noise level
Maxrot=3*pi/18; % maximum head rotation angle
Sim_mode=0; % 0 for ex vivo random rotation, 1 for in vivo limitated rotation
Aceof=-0.02;   
H=zeros(Ndir,3);
FiberDir=zeros([size(mask),Ndir]);    %
if Sim_mode==0
    outpath=fullfile(fpath,strcat('SNR',num2str(SNR),'_Orinum',num2str(Ndir),'_RandomRot'));
    mkdir(outpath);
    save_nii(make_nii(mask),fullfile(outpath,'mask.nii.gz'));
    save_nii(make_nii(WM_mask),fullfile(outpath,'WM_mask.nii.gz'));
    save_nii(make_nii(FAnorm),fullfile(outpath,'FAnorm.nii.gz'));
    for ori=1:Ndir
         theta_x=rand(1)*2*pi-pi;
         theta_y=rand(1)*2*pi-pi;
         Mrot=[1,0,0;0,cos(theta_x),-sin(theta_x);0,sin(theta_x),cos(theta_x)]*[cos(theta_y),0,sin(theta_y);0,1,0;-sin(theta_y),0,cos(theta_y)];
         temp=Mrot*[0,0,1]';
         H(ori,:)=temp';
        temp2=temp(1)*dtipev(:,:,:,1)+temp(2)*dtipev(:,:,:,2)+temp(3)*dtipev(:,:,:,3);   %fiber directions relative to B0 (cos (theta))
        FiberDir(:,:,:,ori)=1-temp2.^2;   %convert to sin(theta)^2.
    end
else
    outpath=fullfile(fpath,strcat('SNR',num2str(SNR),'_Orinum',num2str(Ndir),'_Maxrot',num2str(Maxrot/pi*180)));
    mkdir(outpath);
    save_nii(make_nii(mask),fullfile(outpath,'mask.nii.gz'));
    save_nii(make_nii(WM_mask),fullfile(outpath,'WM_mask.nii.gz'));
    save_nii(make_nii(FAnorm),fullfile(outpath,'FAnorm.nii.gz'));
    H(1,:)=[0,0,1];
    FiberDir(:,:,:,1)=1-dtipev(:,:,:,3).^2;
    for ori=2:Ndir
        theta_r=(2*pi/(Ndir-1))*(ori-2);
        Mrot=[cos(theta_r),-sin(theta_r),0; sin(theta_r),cos(theta_r),0;0,0,1]*[1,0,0;0,cos(Maxrot),-sin(Maxrot);0,sin(Maxrot),cos(Maxrot)];
        temp=Mrot*[0,0,1]';
        H(ori,:)=temp';
        temp2=temp(1)*dtipev(:,:,:,1)+temp(2)*dtipev(:,:,:,2)+temp(3)*dtipev(:,:,:,3); %fiber directions relative to B0 (cos (theta))
        FiberDir(:,:,:,ori)=1-temp2.^2;   %convert to sin(theta)^2.       
    end
end
FiberDir=FiberDir.*WM_mask;
perturb=STI_forward(sim_stitsr,H);
f_non=Aceof*(FiberDir-2/3).*FAnorm.*WM_mask-0.01*WM_mask;    %non-susceptibility-induced related frequency shifts, Eq. [10]
A_GT=Aceof*FAnorm.*WM_mask;
b_GT=-2/3*Aceof*FAnorm.*WM_mask-0.01*WM_mask;  
phi=zeros(size(perturb));
for i=1:size(phi,4)
    temp=perturb(:,:,:,i);
    phi(:,:,:,i)=add_noise(temp,SNR,mask)+f_non(:,:,:,i);
end

%% save simulation data
save_nii(make_nii(phi),fullfile(outpath,'phi.nii.gz'));
save_nii(make_nii(perturb.*mask+f_non),fullfile(outpath,'Cphi.nii.gz'));  %noise-free frequency
save_nii(make_nii(f_non),fullfile(outpath,'f_non_GT.nii.gz'));  %non-susceptibility effect
save_nii(make_nii(FiberDir),fullfile(outpath,'FiberDir.nii.gz'));
save_nii(make_nii(A_GT),fullfile(outpath,'A_GT.nii.gz'));     %microstructure-related coefficient A
save_nii(make_nii(b_GT),fullfile(outpath,'b_GT.nii.gz'));  %orientation-independent part of non-susceptibility-related frequency
save(fullfile(outpath,'H.mat'),'H'); 