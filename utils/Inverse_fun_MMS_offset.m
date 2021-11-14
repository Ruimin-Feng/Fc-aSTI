function [chi_MMS, f_non, flag, relres, iter, resvec] = Inverse_fun_MMS_offset(Params, maxit, tol)
%% inverse function with total variation regularization
%input:
%Params: including the following paramaters
%BrainMask: S1*S2*S3, bainmask;
%phase_tissue: S1*S2*S3*N, normalized field
%OriNum: number of head orientation
%sizeVol: 1*3, matrix size
%H_Matrix: N*0, unit direction vector of main magnetic field in subject frame of reference
%DtiPev: S1*S2*S3*3 principal eigenvector from DTI

%maxit: maximum iteration number
%tol: tolerance

%output:
%chi_par: parallel susceptibility component
%chi_per: perpendicular susceptibility component
%chi_iso: isotropic susceptibility component 
%f_non: non-susceptibility-induced frequency

%% 
OriNum = Params.OriNum;
sizeVol = Params.sizeVol;
VoxNum = prod(sizeVol);     %total voxel number
deltaBArray = Params.phase_tissue;
BrainMask=Params.BrainMask;
Darray=Params.Darray;
I=Params.I;
T1=Params.T1;
MSA=Params.MSA;
%alpha=Params.alpha;
%TV=TVOP;
% calculate MSA-induced phase
MSA_phi=zeros([sizeVol,OriNum]);
 for ori = 1:OriNum
           D_kernel=Darray(:,:,:,:,ori); 
           x_MSA=MSA.*T1.*BrainMask;
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
%b((OriNum*VoxNum+1):end) = 0; 
clear temp deltaBArray

tic
[chi_tensor, flag, relres, iter, resvec] = lsqr(@afun,b,tol,maxit);
toc

% ------------------------------------------------------------------
% change solution format
Chi=reshape(chi_tensor,[sizeVol,2]);
chi_MMS = Chi(:,:,:,1);
f_non = Chi(:,:,:,2);

%resvec = resvec./norm(b);                       % convert to relative residual vector

%% internal function including afun
function res = afun(x,transp_flag)

   if strcmp(transp_flag,'transp')              % y = A'*b
       y = zeros([sizeVol,2]);          % 
       x = reshape(x,[sizeVol,OriNum]);           %input             
       
       for orient_i = 1:OriNum
          Dipole=Darray(:,:,:,1,orient_i)+Darray(:,:,:,5,orient_i)+Darray(:,:,:,9,orient_i); 
          xn=x(:,:,:,orient_i).*BrainMask;              
           y(:,:,:,1)=y(:,:,:,1)+ real(ifftn(fftn(xn).*Dipole)).*BrainMask;
            y(:,:,:,2)=y(:,:,:,2)+xn.*BrainMask;        %
       end
%        temp = alpha*(TV'*(BrainMask.*x(:,:,:,OriNum+1:OriNum+3)));
%        y(:,:,:,2)=y(:,:,:,2)+temp;
       res=y(:);
   elseif strcmp(transp_flag,'notransp')            % y = A*x
        y = zeros([sizeVol,OriNum]);       
        x=reshape(x,[sizeVol,2]);
        
       for orient_i = 1:OriNum
           Dipole=Darray(:,:,:,1,orient_i)+Darray(:,:,:,5,orient_i)+Darray(:,:,:,9,orient_i); 
           x_MMS=x(:,:,:,1).*BrainMask; %
           x_non=x(:,:,:,2).*BrainMask;  %  
           y(:,:,:,orient_i) =real(ifftn(fftn(x_MMS).*Dipole)+x_non).*BrainMask;   
       end
%        TVreg = alpha*BrainMask.*(TV*(x(:,:,:,2) ));
%         y(:,:,:,OriNum+1)=TVreg(:,:,:,1);
%         y(:,:,:,OriNum+2)=TVreg(:,:,:,2);
%         y(:,:,:,OriNum+3)=TVreg(:,:,:,3);
       res=y(:);
    fprintf('+');            
   end
end
end