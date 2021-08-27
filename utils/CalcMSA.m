function []=CalcMSA(path)
%=======================
% This function calculates PEV maps and MSA weighted PEV maps.
% Inputs: 
% path: directory where the input data are located
%=======================
path0=pwd;
cd(path);
name={'_Fc-aSTI', '_aSTI', '_Fc-STI','_STI'};
for i=1:4
    MSA_data=load_nii(['MSA',name{i},'.nii']);
    MSA=MSA_data.img;
    MSA=(MSA-min(MSA(:)))/(max(MSA(:))-min(MSA(:)));
    
    eigenVec_data=load_nii(['Vec',name{i},'.nii']);
    eigenVec=eigenVec_data.img;
    eigenVec3=eigenVec(:,:,:,:,3);

   cMSA=abs(eigenVec3).*repmat(MSA,[1 1 1 3]);
   save_nii(make_nii(cMSA),['cMSA',name{i},'.nii']);
   save_nii(make_nii(eigenVec3),['PEV',name{i},'.nii']);
  
end
 cd(path0);
end







