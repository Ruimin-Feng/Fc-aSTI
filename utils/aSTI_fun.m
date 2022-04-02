    function [ res, tflag ] = aSTI_fun( in, H_Matrix, kx,ky,kz,k2,N_direction,tflag)
        if strcmp(tflag,'transp')
            im = reshape(in, [size(kx,1), N_direction]);
            Res = zeros([size(kx,1), 9]);
            for n = 1:N_direction
                H_Vec = H_Matrix(n,:);
                kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky + H_Vec(3) * kz) ./ (eps + k2);
                %     kH_over_k2=reshape(kH_over_k2,[len,1]);
                Res(:,1) = Res(:,1) + ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* im(:,n);
                
                Res(:,2) = Res(:,2) + ((H_Vec(1)*H_Vec(2))/3 - H_Vec(2)*kx .* kH_over_k2) .* im(:,n);
                
                Res(:,3) = Res(:,3) + ((H_Vec(1)*H_Vec(3))/3 - H_Vec(3)*kx .* kH_over_k2) .* im(:,n);
                
                Res(:,4) = Res(:,4) + ((H_Vec(2)*H_Vec(1))/3 - H_Vec(1)*ky .* kH_over_k2) .* im(:,n);
                
                Res(:,5) = Res(:,5) + ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* im(:,n);
                
                Res(:,6) = Res(:,6) + ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*ky .* kH_over_k2) .* im(:,n);
                
                Res(:,7) = Res(:,7) + ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*kz .* kH_over_k2) .* im(:,n);
                
                Res(:,8) = Res(:,8) + ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*kz .* kH_over_k2) .* im(:,n);
                
                Res(:,9) = Res(:,9) + ((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* im(:,n);
            end
            res = Res(:);
            
        else
            Fx = reshape(in, [size(kx,1),9]);
            Res = zeros([size(kx,1), N_direction]);
            
            for n = 1:N_direction
                H_Vec = H_Matrix(n,:);
                kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky + H_Vec(3) * kz) ./ (eps + k2);
                %   kH_over_k2=reshape(kH_over_k2,[len,1]);
                Res(:,n) = ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,1) + ... % Fx11
                    ((H_Vec(1)*H_Vec(2))/3 -  H_Vec(2)*kx .* kH_over_k2) .* Fx(:,2) + ...    %  Fx12
                    ((H_Vec(1)*H_Vec(3))/3 -  H_Vec(3)*kx .* kH_over_k2) .* Fx(:,3) + ...    %  Fx13
                    ((H_Vec(2)*H_Vec(1))/3 -  H_Vec(1)*ky .* kH_over_k2) .* Fx(:,4) + ...    %  Fx21
                    ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,5) + ...                 %  Fx22
                    ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*ky .*  kH_over_k2) .* Fx(:,6) + ...    %   Fx23
                    ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*kz .*  kH_over_k2) .* Fx(:,7) + ...     %  Fx31
                    ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*kz .*  kH_over_k2) .* Fx(:,8) + ...     %  Fx32
                    ((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Fx(:,9);                        %  Fx33
            end
            res = Res(:);
            fprintf('+')
        end
        
        
    end