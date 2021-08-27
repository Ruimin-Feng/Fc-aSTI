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
[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1,-N(1)/2:N(1)/2-1,-N(3)/2:N(3)/2-1);      % k-space grid
kx = fftshift(kx);      ky = fftshift(ky);      kz = fftshift(kz);
kx=reshape(kx,[N(1)*N(2)*N(3),1]);
ky=reshape(ky,[N(1)*N(2)*N(3),1]);
kz=reshape(kz,[N(1)*N(2)*N(3),1]);
k2 = kx.^2 + ky.^2 + kz.^2;
N_direction = size(phase_tissue, 4);        % number of directions
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


%% solve STI using symmetric+offset model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % number of LSQR iterations    
model='Fc-STI';
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
save_nii(make_nii(chi_res),fullfile(folder,'Fc-STI.nii'));

%% solve STI using asymmetric model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % number of LSQR iterations    
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

%% solve STI using asymmetric+offset model
lsqr_tol = 0.5e-2;                            % LSQR tolerance
lsqr_iter = 100;                             % number of LSQR iterations    
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
end