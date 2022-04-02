function perturb=STI_forward(stitsr,H_Matrix)
%=====================
% This function implements forward STI model
% Inputs:
% stitsr: N1*N2*N3*6, susceptibility tensor elements
% H_Matrix: N_direction*3, each row represents direction vector of the applied magnetic field
%=====================
N_direction=size(H_Matrix,1);
SS=[size(stitsr,1),size(stitsr,2),size(stitsr,3)];
perturb=zeros([SS,N_direction]);
Kstitsr=zeros(size(stitsr));
for i=1:6
    Kstitsr(:,:,:,i)=fftnc(stitsr(:,:,:,i));
end
% k-space grid
[ky,kx,kz] = meshgrid(-SS(2)/2:SS(2)/2-1,-SS(1)/2:SS(1)/2-1,-SS(3)/2:SS(3)/2-1);      
k2 = kx.^2 + ky.^2 + kz.^2;
%forward model
for n = 1:N_direction
    H_Vec = H_Matrix(n,:);
    H_Vec=H_Vec./norm(H_Vec);
    kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky + H_Vec(3) *kz) ./ (eps + k2);
        
    Kperturb= ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Kstitsr(:,:,:,1) + ...                         %   Fx11 
            (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Kstitsr(:,:,:,2) + ...    %   Fx12
            (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*kz + H_Vec(3)*kx) .* kH_over_k2) .* Kstitsr(:,:,:,3) + ...    %   Fx13
            ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Kstitsr(:,:,:,4) + ...                                    %   Fx22
            (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*kz + H_Vec(3)*ky) .* kH_over_k2) .* Kstitsr(:,:,:,5) + ...    %   Fx23
            ((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Kstitsr(:,:,:,6);                                         %   Fx33
    perturb(:,:,:,n)=real(ifftnc(Kperturb));
end
    
end