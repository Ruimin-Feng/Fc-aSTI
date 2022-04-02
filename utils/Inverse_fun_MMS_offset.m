function [chi_MMS, b_fit, flag, relres, iter, resvec] = Inverse_fun_MMS_offset(Params, maxit, tol)
%====================================
%This function is to solve Eq. [A3] in the paper to obtain MMS ans offset maps
% input:
% Params includes:
% OriNum: number of head orientations
% sizeVol: matrix size
% phase_tissue:  H*W*D*NumDir, multi-orientation phase data
% H_Matrix: NumDir*3, direction vector of the magnetic field
% BrainMask: H*W*D, brain mask
% Offset_mask: H*W*D, mask that defines the regions where frequency offset will be considered.
% Darray: H*W*D*9*NumDir, array stores dipole kernel
% T1: T-1/3*I in Eq. [A3]
% MSA: calculated MSA map
% maxit: maximum iteration number
% tol: tolerance
%=========================================
%% 
OriNum = Params.OriNum;
sizeVol = Params.sizeVol;
VoxNum = prod(sizeVol);     %total voxel number
deltaBArray = Params.phase_tissue;
BrainMask=Params.BrainMask;
Offset_mask=Params.Offset_mask;
Darray=Params.Darray;
T1=Params.T1;
MSA=Params.MSA;


% calculate MSA-induced phase
MSA_phi=zeros([sizeVol,OriNum]);
 for ori = 1:OriNum
           D_kernel=Darray(:,:,:,:,ori); 
           x_MSA=MSA.*T1.*Offset_mask;   
           kx_MSA=zeros(sizeVol);     
           for cpt=1:9
               kx_MSA=kx_MSA+fftn(x_MSA(:,:,:,cpt)).*D_kernel(:,:,:,cpt);
           end    
           MSA_phi(:,:,:,ori) =real(ifftn(kx_MSA)).*BrainMask;   
 end
b = zeros((OriNum)*VoxNum, 1); 
for OriInd = 1:OriNum    
    temp = (deltaBArray(:,:,:,OriInd)-MSA_phi(:,:,:,OriInd)).*BrainMask;  
    b(((OriInd - 1)*VoxNum+1) : OriInd*VoxNum) = temp(:);
end
clear temp deltaBArray

tic
[chi_tensor, flag, relres, iter, resvec] = lsqr(@afun,b,tol,maxit);
toc

% ------------------------------------------------------------------
% change solution format
Chi=reshape(chi_tensor,[sizeVol,2]);
chi_MMS = Chi(:,:,:,1);
b_fit = Chi(:,:,:,2);

%% internal function including afun
function res = afun(x,transp_flag)

   if strcmp(transp_flag,'transp')              % y = A'*b
       y = zeros([sizeVol,2]);           
       x = reshape(x,[sizeVol,OriNum]);           %input             
       
       for orient_i = 1:OriNum
          Dipole=Darray(:,:,:,1,orient_i)+Darray(:,:,:,5,orient_i)+Darray(:,:,:,9,orient_i); 
          xn=x(:,:,:,orient_i).*BrainMask;              
           y(:,:,:,1)=y(:,:,:,1)+ real(ifftn(fftn(xn).*Dipole)).*BrainMask;
            y(:,:,:,2)=y(:,:,:,2)+xn.*Offset_mask;      
       end
       res=y(:);
   elseif strcmp(transp_flag,'notransp')            % y = A*x
        y = zeros([sizeVol,OriNum]);       
        x=reshape(x,[sizeVol,2]);    % input
        
       for orient_i = 1:OriNum
           Dipole=Darray(:,:,:,1,orient_i)+Darray(:,:,:,5,orient_i)+Darray(:,:,:,9,orient_i); 
           x_MMS=x(:,:,:,1); %
           x_non=x(:,:,:,2);  %  
           y(:,:,:,orient_i) =real(ifftn(fftn(x_MMS).*Dipole)).*BrainMask++x_non.*Offset_mask;   
       end
       res=y(:);
    fprintf('+');            
   end
end
end