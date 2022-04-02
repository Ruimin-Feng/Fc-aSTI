function [Chi_vec,Chi_eig,MMS,MSA]=Cal_eig(chi_res,msk,model)
%========================
% This function calculates eigenvalue, eigenvector, MMS, and MSA from the
% susceptibility tensor elements
% Inputs:
% chi_res: N1*N2*N3*6 (for symmetric models: STI and Fc-STI) or N1*N2*N3*9
% (for asymmetric models: aSTI and Fc-aSTI)
% msk: N1*N2*N3, tissue mask
% model: 'symm' or 'asymm'
%========================

SS=[size(chi_res,1),size(chi_res,2),size(chi_res,3)];
%modified by Ruimin Feng
Chi_tensor = zeros([SS, 3, 3]);
if strcmp(model,'symm')
    Chi_tensor(:,:,:,1,1) = chi_res(:,:,:,1);
    Chi_tensor(:,:,:,1,2) = chi_res(:,:,:,2);
    Chi_tensor(:,:,:,2,1) = chi_res(:,:,:,2);
    Chi_tensor(:,:,:,1,3) = chi_res(:,:,:,3);
    Chi_tensor(:,:,:,3,1) = chi_res(:,:,:,3);
    Chi_tensor(:,:,:,2,2) = chi_res(:,:,:,4);  %
    Chi_tensor(:,:,:,2,3) = chi_res(:,:,:,5);   %
    Chi_tensor(:,:,:,3,2) = chi_res(:,:,:,5);   %
    Chi_tensor(:,:,:,3,3) = chi_res(:,:,:,6);   %
else
    Chi_tensor(:,:,:,1,1) = chi_res(:,:,:,1);
    Chi_tensor(:,:,:,1,2) = chi_res(:,:,:,2);
    Chi_tensor(:,:,:,2,1) = chi_res(:,:,:,4);
    Chi_tensor(:,:,:,1,3) = chi_res(:,:,:,3);
    Chi_tensor(:,:,:,3,1) = chi_res(:,:,:,7);
    Chi_tensor(:,:,:,2,2) = chi_res(:,:,:,5);  %
    Chi_tensor(:,:,:,2,3) = chi_res(:,:,:,6);   %
    Chi_tensor(:,:,:,3,2) = chi_res(:,:,:,8);   %
    Chi_tensor(:,:,:,3,3) = chi_res(:,:,:,9);   
end
    

mask_tensor = msk(:);

Chi_tensor = permute(Chi_tensor, [4,5,1,2,3]);
chi_tensor = reshape(Chi_tensor, [3,3, numel(mask_tensor)]);

Chi_eig = zeros(numel(mask_tensor),3);
Chi_vec=zeros(numel(mask_tensor),3,3);

tic
for v = 1:length(mask_tensor)
    if mask_tensor(v) ~= 0
        [V,D] = eig(chi_tensor(:,:,v));
        Chi_eig(v,:) = diag(D)';
        Chi_vec(v,:,:)=V;
    end
end
toc

Chi_eig = reshape(Chi_eig, [SS, 3]);
Chi_vec=reshape(Chi_vec,[SS,3,3]);
MMS = mean(Chi_eig,4);
MSA = Chi_eig(:,:,:,3) - (Chi_eig(:,:,:,1) + Chi_eig(:,:,:,2)) / 2; 
end