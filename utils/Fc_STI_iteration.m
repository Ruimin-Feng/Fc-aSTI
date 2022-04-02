function [chi_MMS,chi_MSA,Pev,Offset]= Fc_STI_iteration(tissue_phase,H_Matrix,mask,offset_mask,iter_num,outpath)
%====================================
%This function is to solve Fc_STI  according to the four steps of Figure 1 in the paper:
% An improved asymmetric susceptibility tensor imaging model with frequency offset correction
% (Feng et al)
%tissue_phase: H*W*D*NumDir, multi-orientation phase data
%H_Matrix: NumDir*3, direction vector of the magnetic field
%mask: H*W*D, brain mask
%offset_mask: H*W*D, mask that defines the regions where frequency offset will be considered.
%iter_num: number of iterations, 0 for no iteration
%outpath: outpath directory
%=========================================
%% 
N_direction = size(tissue_phase,4);    %number of head direction
sizeVol =size(mask);                           %volume size
%% parameter preparation
[KY_Grid, KX_Grid, KZ_Grid] = meshgrid(-sizeVol(2)/2:sizeVol(2)/2-1,-sizeVol(1)/2:sizeVol(1)/2-1,-sizeVol(3)/2:sizeVol(3)/2-1);      % k-space grid
KX_Grid = fftshift(KX_Grid);      KY_Grid = fftshift(KY_Grid);      KZ_Grid = fftshift(KZ_Grid);
KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;          % k^2
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
Darray=zeros([sizeVol,9,N_direction]);     %dipole kernel
for orient_i = 1:N_direction
    H_Vec = H_Matrix(orient_i,:);
    HTH_Vec=HTH(orient_i,:);
    HTH_Mat=repmat(reshape(HTH_Vec,1,1,1,9),[sizeVol,1]);
    kH_over_k2 = (H_Vec(1) * KX_Grid + H_Vec(2) * KY_Grid + H_Vec(3) * KZ_Grid) ./ (eps + KSq);
    %KTH
    KTH=zeros([sizeVol,9]);
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

%% step 1: Solving seven unknowns using lsqr.
lsqr_tol = 0.001;                            % LSQR tolerance
lsqr_iter = 30;
[lsqr_out, ~, ~, ~] = lsqr(@afun, tissue_phase(:), lsqr_tol, lsqr_iter);
out=reshape(lsqr_out,[sizeVol,7]);

%% step 2: Estimate MSA and PEV
chi_tensor=out(:,:,:,1:6);
[Chi_vec,~,~,chi_MSA]=Cal_eig(chi_tensor,mask,'symm');
Pev=Chi_vec(:,:,:,:,3);
save_nii(make_nii(chi_MSA),fullfile(outpath,'MSA_Fc-STI_iter0.nii'));
save_nii(make_nii(abs(Pev)),fullfile(outpath,'absPEV_Fc-STI_iter0.nii')); 
save_nii(make_nii(Pev),fullfile(outpath,'PEV_Fc-STI_iter0.nii'));   

%% step 3: Solving MMS and offset map based on the cylindrical symmetry constraint imposed on the tensor
maxit=30;
tol=1e-3;
T1=zeros([sizeVol,9]);
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
Params.MSA=chi_MSA;        
Params.OriNum=N_direction;
Params.sizeVol=sizeVol;
Params.phase_tissue=tissue_phase;
Params.BrainMask=mask;
Params.Darray=Darray;
Params.Offset_mask=offset_mask;
[chi_MMS, Offset, flag, relres, iter, ~] = Inverse_fun_MMS_offset(Params, maxit, tol);
disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)]);   
save_nii(make_nii(chi_MMS),fullfile(outpath,'MMS_Fc-STI_iter0.nii'));
save_nii(make_nii(Offset),fullfile(outpath,'Fres_Fc-STI_iter0.nii'));

%% Step 4: iteration process
if iter_num>0
    for cpt=1:iter_num
        % solve MSA and PEV
        phase_temp=(tissue_phase-Offset).*mask;
        [~,~,chi_MSA,Pev]=STI_inverse(phase_temp,H_Matrix,mask);
        save_nii(make_nii(chi_MSA),fullfile(outpath,['MSA_Fc-STI_iter',num2str(cpt),'.nii']));
        save_nii(make_nii(Pev),fullfile(outpath,['PEV_Fc-STI_iter',num2str(cpt),'.nii']));
        save_nii(make_nii(abs(Pev)),fullfile(outpath,['absPEV_Fc-STI_iter',num2str(cpt),'.nii']));
        % solve MMS and Offset
        maxit=30;
        tol=1e-3;
        T1=zeros([sizeVol,9]);
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
        Params.MSA=chi_MSA;
        Params.OriNum=N_direction;
        Params.sizeVol=sizeVol;
        Params.phase_tissue=tissue_phase;
        Params.BrainMask=mask;
        Params.Darray=Darray;
        Params.Offset_mask=offset_mask;
        [chi_MMS, Offset, flag, relres, iter, ~] = Inverse_fun_MMS_offset(Params, maxit, tol);
        disp(['Flag: ', num2str(flag), '  Relres: ', num2str(relres), '  Iter: ', num2str(iter)]);
         save_nii(make_nii(chi_MMS),fullfile(outpath,['MMS_Fc-STI_iter',num2str(cpt),'.nii']));
        save_nii(make_nii(Offset),fullfile(outpath,['Offset_Fc-STI_iter',num2str(cpt),'.nii']));    
    end
end





    function res = afun( in,tflag)
        if strcmp(tflag,'transp')
            im = reshape(in, [sizeVol, N_direction]);
            Res = zeros([sizeVol, 7]);
            for n = 1:N_direction
                Res(:,:,:,1) = Res(:,:,:,1) + real(ifftn(Darray(:,:,:,1,n) .*fftn( im(:,:,:,n).*mask)));
                
                Res(:,:,:,2) = Res(:,:,:,2) +real(ifftn( (Darray(:,:,:,2,n)+Darray(:,:,:,4,n)) .*fftn( im(:,:,:,n).*mask)));
                
                Res(:,:,:,3) = Res(:,:,:,3) + real(ifftn((Darray(:,:,:,3,n)+Darray(:,:,:,7,n)) .*fftn( im(:,:,:,n).*mask)));
                
                Res(:,:,:,4) = Res(:,:,:,4) +real(ifftn( Darray(:,:,:,5,n) .* fftn( im(:,:,:,n).*mask)));
                
                Res(:,:,:,5) = Res(:,:,:,5) +real(ifftn((Darray(:,:,:,6,n)+Darray(:,:,:,8,n)) .*fftn( im(:,:,:,n).*mask)));
                
                Res(:,:,:,6) = Res(:,:,:,6) +real(ifftn( Darray(:,:,:,9,n).* fftn( im(:,:,:,n).*mask)));
              
                Res(:,:,:,7)=Res(:,:,:,7)+offset_mask.*im(:,:,:,n);
            end
            res = Res(:);
            
        else
            im = reshape(in, [sizeVol,7]);
            Res = zeros([sizeVol, N_direction]);
            
            for n = 1:N_direction               
                Res(:,:,:,n) =mask.*real(ifftn( Darray(:,:,:,1,n) .*fftn( im(:,:,:,1)) + ...
                    (Darray(:,:,:,2,n)+Darray(:,:,:,4,n)) .* fftn(im(:,:,:,2)) + ...    
                    (Darray(:,:,:,3,n)+Darray(:,:,:,7,n)) .*fftn( im(:,:,:,3)) + ...   
                    Darray(:,:,:,5,n) .* fftn(im(:,:,:,4)) + ...  
                    (Darray(:,:,:,6,n)+Darray(:,:,:,8,n)) .*fftn( im(:,:,:,5)) + ...                
                    Darray(:,:,:,9,n) .*fftn( im(:,:,:,6)))) + ...   
                    offset_mask.*im(:,:,:,7);                                                                                                       %  Offset
            end
            res = Res(:);
            fprintf('+')
        end
    end
end

