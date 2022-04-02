function [chi_tensor,chi_MMS,chi_MSA,Pev]=aSTI_inverse(phase_tissue,H_Matrix,mask,varargin)       
%====================================
% This function is to solve the aSTI 
% Parallel Computing Toolbox is used to speed up inversion process.
% phase_tissue: H*W*D*NumDir, multi-orientation phase data
% H_Matrix: NumDir*3, direction vector of the magnetic field
% mask: H*W*D, brain mask
%=========================================N=size(mask);
%% prepare parameters
N=size(mask);
N_direction = size(phase_tissue, 4);        % no of directions
[KY_Grid, KX_Grid, KZ_Grid] = meshgrid(-N(2)/2:N(2)/2-1,-N(1)/2:N(1)/2-1,-N(3)/2:N(3)/2-1);      % k-space grid
KX_Grid = fftshift(KX_Grid);      KY_Grid = fftshift(KY_Grid);      KZ_Grid = fftshift(KZ_Grid);
kx=reshape(KX_Grid,[N(1)*N(2)*N(3),1]);
ky=reshape(KY_Grid,[N(1)*N(2)*N(3),1]);
kz=reshape(KZ_Grid,[N(1)*N(2)*N(3),1]);
k2 = kx.^2 + ky.^2 + kz.^2;
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

%%
lsqr_tol = 0.005;                            % LSQR tolerance
lsqr_iter = 100;                             % no of LSQR iterations      %*************
tic
parfor loopvar = 1:poolnum
    warning off all
    [respar{loopvar}, flag, relres, iter] = lsqr(@aSTI_fun, thetai0par{loopvar}, lsqr_tol, lsqr_iter, [], [], [], H_Matrix,kxpar{loopvar},kypar{loopvar},kzpar{loopvar},k2par{loopvar},N_direction);  
    disp(['Solving aSTI---Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)])
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

chi_tensor = zeros(size(Fchi_res));
for ind = 1:9
    chi_tensor(:,:,:,ind) = real( ifftn(Fchi_res(:,:,:,ind)) ) .* mask;             % susceptibility tensor in image space
end
[symmpart,~]=tensor_decomp(chi_tensor);
[Chi_vec,~,chi_MMS,chi_MSA]=Cal_eig(symmpart,mask,'asymm');
Pev=Chi_vec(:,:,:,:,3);
if nargin>3
    try
        save_nii(make_nii(chi_MMS),fullfile(varargin{1},'MMS_aSTI.nii'));
        save_nii(make_nii(chi_MSA),fullfile(varargin{1},'MSA_aSTI.nii'));
        save_nii(make_nii(abs(Pev)),fullfile(varargin{1},'absPEV_aSTI.nii')); 
        save_nii(make_nii(Pev),fullfile(varargin{1},'PEV_aSTI.nii'));    
    catch
        disp('Failed to save the results');
    end
end
    

end