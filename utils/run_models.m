function []=run_models(phase_tissue, H_Matrix, mask, folder)
%============================
% This function runs STI, Fc-STI, aSTI, and Fc-aSTI models. 
% Parallel Computing Toolbox is used to speed up inversion process.
% Inputs: 
% phase_tissue: N1*N2*N3*N_direction, phase data
% H_Matrix: N_direction*3, each row represents unit direction vector of the applied magnetic field
% mask: N1*N2*N3, tissue mask
% folder: output directory
%=============================
%% prepare parameters
N=size(mask);
[KY_Grid, KX_Grid, KZ_Grid] = meshgrid(-N(2)/2:N(2)/2-1,-N(1)/2:N(1)/2-1,-N(3)/2:N(3)/2-1);      % k-space grid
KX_Grid = fftshift(KX_Grid);      KY_Grid = fftshift(KY_Grid);      KZ_Grid = fftshift(KZ_Grid);
kx=reshape(KX_Grid,[N(1)*N(2)*N(3),1]);
ky=reshape(KY_Grid,[N(1)*N(2)*N(3),1]);
kz=reshape(KZ_Grid,[N(1)*N(2)*N(3),1]);
k2 = kx.^2 + ky.^2 + kz.^2;
N_direction = size(phase_tissue, 4);        % no of directions
Phase_tissue = zeros(size(phase_tissue));
for ind = 1:N_direction
    Phase_tissue(:,:,:,ind) = fftn(phase_tissue(:,:,:,ind));
end
Phase_tissue=reshape(Phase_tissue,[N(1)*N(2)*N(3),N_direction]);
poolnum=12;
sizepar=round(N(1)*N(2)*N(3)/poolnum);
kxpar=cell(poolnum,1);
kypar=cell(poolnum,1);
kzpar=cell(poolnum,1);
k2par=cell(poolnum,1);
thetai0par=cell(poolnum,1);
for i=1:poolnum
    startpar=(i-1)*sizepar+1;
    if i==poolnum
        endpar=N(1)*N(2)*N(3);
    else
        endpar=i*sizepar;
    end
    kxpar{i}=kx(startpar:endpar,:);
    kypar{i}=ky(startpar:endpar,:);
    kzpar{i}=kz(startpar:endpar,:);
    k2par{i}=k2(startpar:endpar,:);
    temp=Phase_tissue(startpar:endpar,:);
    thetai0par{i}=temp(:);
end
disp('Parallel variables created')

%% prepare parameters for solving second step
KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;  
HTH=zeros([N_direction,9]);
HTH(:,1)=H_Matrix(:,1).*H_Matrix(:,1);
HTH(:,2)=H_Matrix(:,1).*H_Matrix(:,2);
HTH(:,3)=H_Matrix(:,1).*H_Matrix(:,3);
HTH(:,4)=H_Matrix(:,2).*H_Matrix(:,1);
HTH(:,5)=H_Matrix(:,2).*H_Matrix(:,2);
HTH(:,6)=H_Matrix(:,2).*H_Matrix(:,3);
HTH(:,7)=H_Matrix(:,3).*H_Matrix(:,1);
HTH(:,8)=H_Matrix(:,3).*H_Matrix(:,2);
HTH(:,9)=H_Matrix(:,3).*H_Matrix(:,3);
Darray=zeros([N,9,N_direction]);     %dipole kernel
for orient_i = 1:N_direction
    H_Vec = H_Matrix(orient_i,:);
    HTH_Vec=HTH(orient_i,:);
    HTH_Mat=repmat(reshape(HTH_Vec,1,1,1,9),[N,1]);
    kH_over_k2 = (H_Vec(1) * KX_Grid + H_Vec(2) * KY_Grid + H_Vec(3) * KZ_Grid) ./ (eps + KSq);
    %KTH
    KTH=zeros([N,9]);
    KTH(:,:,:,1)=KX_Grid.*H_Vec(1);
    KTH(:,:,:,2)=KX_Grid.*H_Vec(2);
    KTH(:,:,:,3)=KX_Grid.*H_Vec(3);
    KTH(:,:,:,4)=KY_Grid.*H_Vec(1);
    KTH(:,:,:,5)=KY_Grid.*H_Vec(2);
    KTH(:,:,:,6)=KY_Grid.*H_Vec(3);
    KTH(:,:,:,7)=KZ_Grid.*H_Vec(1);
    KTH(:,:,:,8)=KZ_Grid.*H_Vec(2);
    KTH(:,:,:,9)=KZ_Grid.*H_Vec(3);
    Darray(:,:,:,:,orient_i)=1/3*HTH_Mat-kH_over_k2.*KTH;
end
Params.OriNum = N_direction;
Params.sizeVol = N;
Params.H_Matrix = H_Matrix;
Params.phase_tissue = phase_tissue;
Params.BrainMask=mask;
Params.Darray=Darray;

%% solve STI using symmetric model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % number of LSQR iterations    
model='STI';
tic
parfor loopvar = 1:poolnum
    warning off all
    [respar{loopvar}, flag, relres, iter] = lsqr(@compute_STI_fast, thetai0par{loopvar}, lsqr_tol, lsqr_iter, [], [], [],H_Matrix, kxpar{loopvar},kypar{loopvar},kzpar{loopvar},k2par{loopvar},N,N_direction,model);  
    disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
end
toc

Fchi_res=zeros([N(1)*N(2)*N(3),6]);
for i=1:poolnum
    startpar=(i-1)*sizepar+1;
    if i==poolnum
        endpar=N(1)*N(2)*N(3);
    else
        endpar=i*sizepar;
    end
    Fchi_res(startpar:endpar,:)=reshape(respar{i},[],6);
end
Fchi_res = reshape(Fchi_res, [N,6]);             % susceptibility tensor in k-space

chi_res = zeros(size(Fchi_res));
for ind = 1:6
    chi_res(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
save_nii(make_nii(chi_res),fullfile(folder,'STI.nii'));
[Chi_vec,~,MMS,MSA]=Cal_eig(chi_res,mask,'symm');
save_nii(make_nii(MMS),fullfile(folder,'MMS_STI.nii'));
save_nii(make_nii(MSA),fullfile(folder,'MSA_STI.nii'));
save_nii(make_nii(Chi_vec(:,:,:,:,3)),fullfile(folder,'PEV_STI.nii'));
save_nii(make_nii(abs(Chi_vec(:,:,:,:,3))),fullfile(folder,'absPEV_STI.nii'));

%% solve STI using symmetric+offset model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % no of LSQR iterations    
model='Fc-STI';
maxit=100;
tol=1e-3;
tic
parfor loopvar = 1:poolnum
    warning off all
    [respar{loopvar}, flag, relres, iter] = lsqr(@compute_STI_fast, thetai0par{loopvar}, lsqr_tol, lsqr_iter, [], [], [],H_Matrix, kxpar{loopvar},kypar{loopvar},kzpar{loopvar},k2par{loopvar},N,N_direction,model);  
    disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
end
toc

Fchi_res=zeros([N(1)*N(2)*N(3),7]);
for i=1:poolnum
    startpar=(i-1)*sizepar+1;
    if i==poolnum
        endpar=N(1)*N(2)*N(3);
    else
        endpar=i*sizepar;
    end
    Fchi_res(startpar:endpar,:)=reshape(respar{i},[],7);
end
Fchi_res = reshape(Fchi_res, [N,7]);             % susceptibility tensor in k-space

chi_res = zeros(size(Fchi_res));
for ind = 1:7
    chi_res(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
[Chi_vec,~,~,MSA]=Cal_eig(chi_res(:,:,:,1:6),mask,'symm');
save_nii(make_nii(MSA),fullfile(folder,'MSA_Fc-STI.nii'));
save_nii(make_nii(Chi_vec(:,:,:,:,3)),fullfile(folder,'PEV_Fc-STI.nii'));
save_nii(make_nii(abs(Chi_vec(:,:,:,:,3))),fullfile(folder,'absPEV_Fc-STI.nii'));
%second step
Pev=Chi_vec(:,:,:,:,3);
T1=zeros([N,9]);
I=zeros([N,9]);
V1=Pev(:,:,:,1);
V2=Pev(:,:,:,2);
V3=Pev(:,:,:,3);
T1(:,:,:,1)=V1.*V1-1/3;
T1(:,:,:,2)=V1.*V2;
T1(:,:,:,3)=V1.*V3;
T1(:,:,:,4)=V2.*V1;
T1(:,:,:,5)=V2.*V2-1/3;
T1(:,:,:,6)=V2.*V3;
T1(:,:,:,7)=V3.*V1;
T1(:,:,:,8)=V3.*V2;
T1(:,:,:,9)=V3.*V3-1/3;
Params.T1=T1.*mask;   

I(:,:,:,1)=1;
I(:,:,:,2)=0;
I(:,:,:,3)=0;
I(:,:,:,4)=0;
I(:,:,:,5)=1;
I(:,:,:,6)=0;
I(:,:,:,7)=0;
I(:,:,:,8)=0; 
I(:,:,:,9)=1;
Params.I=I.*mask;  
Params.MSA=MSA;
[chi_MMS, f_non, flag, relres, iter] = Inverse_fun_MMS_offset(Params, maxit, tol);
disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)]);
save_nii(make_nii(chi_MMS),fullfile(folder,'MMS_Fc-STI.nii'));
save_nii(make_nii(f_non),fullfile(folder,'f_non_Fc-STI.nii'));

%% solve STI using asymmetric model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % no of LSQR iterations    
model='aSTI';
tic
parfor loopvar = 1:poolnum
    warning off all
    [respar{loopvar}, flag, relres, iter] = lsqr(@compute_STI_fast, thetai0par{loopvar}, lsqr_tol, lsqr_iter, [], [], [],H_Matrix, kxpar{loopvar},kypar{loopvar},kzpar{loopvar},k2par{loopvar},N,N_direction,model);  
    disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
end
toc

Fchi_res=zeros([N(1)*N(2)*N(3),9]);
for i=1:poolnum
    startpar=(i-1)*sizepar+1;
    if i==poolnum
        endpar=N(1)*N(2)*N(3);
    else
        endpar=i*sizepar;
    end
    Fchi_res(startpar:endpar,:)=reshape(respar{i},[],9);
end
Fchi_res = reshape(Fchi_res, [N,9]);             % susceptibility tensor in k-space

chi_res = zeros(size(Fchi_res));
for ind = 1:9
    chi_res(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
save_nii(make_nii(chi_res),fullfile(folder,'aSTI.nii'));
[asymm_symmpart,asymm_antisymmpart]=tensor_decomp(chi_res);
save_nii(make_nii(asymm_symmpart),fullfile(folder,'aSTI_symmpart.nii'));
save_nii(make_nii(asymm_antisymmpart),fullfile(folder,'aSTI_antisymmpart.nii'));
[Chi_vec,~,MMS,MSA]=Cal_eig(asymm_symmpart,mask,'asymm');
save_nii(make_nii(MMS),fullfile(folder,'MMS_aSTI.nii'));
save_nii(make_nii(MSA),fullfile(folder,'MSA_aSTI.nii'));
save_nii(make_nii(Chi_vec(:,:,:,:,3)),fullfile(folder,'PEV_aSTI.nii'));
save_nii(make_nii(abs(Chi_vec(:,:,:,:,3))),fullfile(folder,'absPEV_aSTI.nii'));

%% solve STI using asymmetric+offset model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % no of LSQR iterations    
model='Fc-aSTI';
tic
parfor loopvar = 1:poolnum
    warning off all
    [respar{loopvar}, flag, relres, iter] = lsqr(@compute_STI_fast, thetai0par{loopvar}, lsqr_tol, lsqr_iter, [], [], [],H_Matrix, kxpar{loopvar},kypar{loopvar},kzpar{loopvar},k2par{loopvar},N,N_direction,model);  
    disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
end
toc

Fchi_res=zeros([N(1)*N(2)*N(3),10]);
for i=1:poolnum
    startpar=(i-1)*sizepar+1;
    if i==poolnum
        endpar=N(1)*N(2)*N(3);
    else
        endpar=i*sizepar;
    end
    Fchi_res(startpar:endpar,:)=reshape(respar{i},[],10);
end
Fchi_res = reshape(Fchi_res, [N,10]);             % susceptibility tensor in k-space

chi_res = zeros(size(Fchi_res));
for ind = 1:10
    chi_res(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
save_nii(make_nii(chi_res),fullfile(folder,'Fc-aSTI.nii'));
[asymm_Fs_symmpart,asymm_Fs_antisymmpart]=tensor_decomp(chi_res(:,:,:,1:9));
save_nii(make_nii(asymm_Fs_symmpart),fullfile(folder,'Fc-aSTI_symmpart.nii'));
save_nii(make_nii(asymm_Fs_antisymmpart),fullfile(folder,'Fc-aSTI_antisymmpart.nii'));
[Chi_vec,~,~,MSA]=Cal_eig(asymm_Fs_symmpart,mask,'asymm');
save_nii(make_nii(MSA),fullfile(folder,'MSA_Fc-aSTI.nii'));
save_nii(make_nii(Chi_vec(:,:,:,:,3)),fullfile(folder,'PEV_Fc-aSTI.nii'));
save_nii(make_nii(abs(Chi_vec(:,:,:,:,3))),fullfile(folder,'absPEV_Fc-aSTI.nii'));

%second step
maxit=100;
tol=1e-3;
Pev=Chi_vec(:,:,:,:,3);
T1=zeros([N,9]);
I=zeros([N,9]);
V1=Pev(:,:,:,1);
V2=Pev(:,:,:,2);
V3=Pev(:,:,:,3);
T1(:,:,:,1)=V1.*V1-1/3;
T1(:,:,:,2)=V1.*V2;
T1(:,:,:,3)=V1.*V3;
T1(:,:,:,4)=V2.*V1;
T1(:,:,:,5)=V2.*V2-1/3;
T1(:,:,:,6)=V2.*V3;
T1(:,:,:,7)=V3.*V1;
T1(:,:,:,8)=V3.*V2;
T1(:,:,:,9)=V3.*V3-1/3;
Params.T1=T1.*mask;   

I(:,:,:,1)=1;
I(:,:,:,2)=0;
I(:,:,:,3)=0;
I(:,:,:,4)=0;
I(:,:,:,5)=1;
I(:,:,:,6)=0;
I(:,:,:,7)=0;
I(:,:,:,8)=0; 
I(:,:,:,9)=1;
Params.I=I.*mask;  
Params.MSA=MSA;
[chi_MMS, f_non, flag, relres, iter] = Inverse_fun_MMS_offset(Params, maxit, tol);
disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)]);
save_nii(make_nii(chi_MMS),fullfile(folder,'MMS_Fc-aSTI.nii'));
save_nii(make_nii(f_non),fullfile(folder,'f_non_Fc-aSTI.nii'));
end