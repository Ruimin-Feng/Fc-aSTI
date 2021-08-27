function [ res, tflag ] = compute_STI_fast( in, H_Matrix, kx,ky,kz,k2,N,N_direction,model, tflag)
 %===========================
 % This function implements STI, Fc-STI, aSTI, and Fc-aSTI reconstruction
 % Inputs: 
 % in: phase data
 % H_Matrix: direction vector of the magnetic field
 % kx,ky,kz: Fourier domain coordinates 
 % k2: k2=kx^2+ky^2+kz^2
 % N: size of phase data
 % N_direction: number of head orientations
 % model: reconstruction model(STI, Fc-STI, aSTI, or Fc-aSTI) 
 %===========================
 
 params = [];
 params.SS = N;
 params.N_direction = N_direction;
 params.H_Matrix = H_Matrix;
 params.kx=kx;
 params.ky=ky;
 params.kz=kz;
 params.k2 = k2;
 params.model=model;
%% STI model
if strcmp(params.model,'STI')
    if strcmp(tflag,'transp')
        im = reshape(in, size(kx,1), params.N_direction);
        Res = zeros([size(kx,1), 6]);
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
        %    kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,1) = Res(:,1) + ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,2) = Res(:,2) + (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*params.ky + H_Vec(2)*params.kx) .* kH_over_k2) .* im(:,n);
            
            Res(:,3) = Res(:,3) + (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*params.kz + H_Vec(3)*params.kx) .* kH_over_k2) .* im(:,n);
            
            Res(:,4) = Res(:,4) + ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,5) = Res(:,5) + (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*params.kz + H_Vec(3)*params.ky) .* kH_over_k2) .* im(:,n);
            
            Res(:,6) = Res(:,6) + ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* im(:,n);
        end
        res = Res(:);
   
    else
        Fx = reshape(in, [size(kx,1),6]);
        Res = zeros([size(kx,1), params.N_direction]);
        
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
        %    kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,n) = ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* Fx(:,1) + ...                         %   Fx11
                (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*params.ky + H_Vec(2)*params.kx) .* kH_over_k2) .* Fx(:,2) + ...    %   Fx12
                (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*params.kz + H_Vec(3)*params.kx) .* kH_over_k2) .* Fx(:,3) + ...    %   Fx13
                ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* Fx(:,4) + ...                                    %   Fx22
                (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*params.kz + H_Vec(3)*params.ky) .* kH_over_k2) .* Fx(:,5) + ...    %   Fx23
                ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* Fx(:,6);                                         %   Fx33
        end   
    res = Res(:);   
    fprintf('+')
    end


%% Fc-STI model  
elseif strcmp(params.model,'Fc-STI')
    if strcmp(tflag,'transp')
        im = reshape(in, [size(kx,1), params.N_direction]);
        Res = zeros([size(kx,1), 7]);
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
      %      kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,1) = Res(:,1) + ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,2) = Res(:,2) + (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*params.ky + H_Vec(2)*params.kx) .* kH_over_k2) .* im(:,n);
            
            Res(:,3) = Res(:,3) + (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*params.kz + H_Vec(3)*params.kx) .* kH_over_k2) .* im(:,n);
            
            Res(:,4) = Res(:,4) + ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,5) = Res(:,5) + (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*params.kz + H_Vec(3)*params.ky) .* kH_over_k2) .* im(:,n);
            
            Res(:,6) = Res(:,6) + ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* im(:,n);
            
            Res(:,7)=Res(:,7)+im(:,n);
        end
        res = Res(:);
   
    else
        Fx = reshape(in, [size(kx,1),7]);
        Res = zeros([size(kx,1), params.N_direction]);
        
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
        %    kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,n) = ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* Fx(:,1) + ...                         %   Fx11
                (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*params.ky + H_Vec(2)*params.kx) .* kH_over_k2) .* Fx(:,2) + ...    %   Fx12
                (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*params.kz + H_Vec(3)*params.kx) .* kH_over_k2) .* Fx(:,3) + ...    %   Fx13
                ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* Fx(:,4) + ...                                    %   Fx22
                (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*params.kz + H_Vec(3)*params.ky) .* kH_over_k2) .* Fx(:,5) + ...    %   Fx23
                ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* Fx(:,6)+Fx(:,7);                                         %   Fx33
        end   
    res = Res(:);   
    fprintf('+')
    end

%% aSTI model
elseif strcmp(params.model,'aSTI')
    if strcmp(tflag,'transp')
        im = reshape(in, [size(kx,1), params.N_direction]);
        Res = zeros([size(kx,1), 9]);
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
       %     kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,1) = Res(:,1) + ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,2) = Res(:,2) + ((H_Vec(1)*H_Vec(2))/3 - H_Vec(2)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,3) = Res(:,3) + ((H_Vec(1)*H_Vec(3))/3 - H_Vec(3)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,4) = Res(:,4) + ((H_Vec(2)*H_Vec(1))/3 - H_Vec(1)*params.ky .* kH_over_k2) .* im(:,n);
                       
            Res(:,5) = Res(:,5) + ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,6) = Res(:,6) + ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,7) = Res(:,7) + ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*params.kz .* kH_over_k2) .* im(:,n);
            
            Res(:,8) = Res(:,8) + ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*params.kz .* kH_over_k2) .* im(:,n);
            
            Res(:,9) = Res(:,9) + ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* im(:,n);
        end
        res = Res(:);
   
    else
        Fx = reshape(in, [size(kx,1),9]);
        Res = zeros([size(kx,1), params.N_direction]);
        
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
         %   kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,n) = ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* Fx(:,1) + ... % Fx11
                ((H_Vec(1)*H_Vec(2))/3 -  H_Vec(2)*params.kx .* kH_over_k2) .* Fx(:,2) + ...    %  Fx12
                ((H_Vec(1)*H_Vec(3))/3 -  H_Vec(3)*params.kx .* kH_over_k2) .* Fx(:,3) + ...    %  Fx13
                ((H_Vec(2)*H_Vec(1))/3 -  H_Vec(1)*params.ky .* kH_over_k2) .* Fx(:,4) + ...    %  Fx21
                ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* Fx(:,5) + ...                 %  Fx22
                ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*params.ky .*  kH_over_k2) .* Fx(:,6) + ...    %   Fx23
                ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*params.kz .*  kH_over_k2) .* Fx(:,7) + ...     %  Fx31
                ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*params.kz .*  kH_over_k2) .* Fx(:,8) + ...     %  Fx32
                ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* Fx(:,9);                        %  Fx33
        end   
    res = Res(:);   
    fprintf('+')
    end
    
    
    
    %% Fc-aSTI model
elseif strcmp(params.model,'Fc-aSTI')
   if strcmp(tflag,'transp')
        im = reshape(in, [size(kx,1), params.N_direction]);
        Res = zeros([size(kx,1), 10]);
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
         %   kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,1) = Res(:,1) + ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,2) = Res(:,2) + ((H_Vec(1)*H_Vec(2))/3 - H_Vec(2)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,3) = Res(:,3) + ((H_Vec(1)*H_Vec(3))/3 - H_Vec(3)*params.kx .* kH_over_k2) .* im(:,n);
            
            Res(:,4) = Res(:,4) + ((H_Vec(2)*H_Vec(1))/3 - H_Vec(1)*params.ky .* kH_over_k2) .* im(:,n);
                       
            Res(:,5) = Res(:,5) + ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,6) = Res(:,6) + ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*params.ky .* kH_over_k2) .* im(:,n);
            
            Res(:,7) = Res(:,7) + ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*params.kz .* kH_over_k2) .* im(:,n);
            
            Res(:,8) = Res(:,8) + ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*params.kz .* kH_over_k2) .* im(:,n);
            
            Res(:,9) = Res(:,9) + ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* im(:,n);
            Res(:,10) = Res(:,10) +im(:,n);
        end
        res = Res(:);
   
    else
        Fx = reshape(in, [size(kx,1),10]);
        Res = zeros([size(kx,1), params.N_direction]);
        
        for n = 1:params.N_direction
            H_Vec = params.H_Matrix(n,:);
            kH_over_k2 = (H_Vec(1) * params.kx + H_Vec(2) * params.ky + H_Vec(3) * params.kz) ./ (eps + params.k2);
          %  kH_over_k2=reshape(kH_over_k2,[len,1]);
            Res(:,n) = ((H_Vec(1)^2)/3 - H_Vec(1)*params.kx .* kH_over_k2) .* Fx(:,1) + ... % Fx11
                ((H_Vec(1)*H_Vec(2))/3 -  H_Vec(2)*params.kx .* kH_over_k2) .* Fx(:,2) + ...    %  Fx12
                ((H_Vec(1)*H_Vec(3))/3 -  H_Vec(3)*params.kx .* kH_over_k2) .* Fx(:,3) + ...    %  Fx13
                ((H_Vec(2)*H_Vec(1))/3 -  H_Vec(1)*params.ky .* kH_over_k2) .* Fx(:,4) + ...    %  Fx21
                ((H_Vec(2)^2)/3 - H_Vec(2)*params.ky .* kH_over_k2) .* Fx(:,5) + ...                 %  Fx22
                ((H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*params.ky .*  kH_over_k2) .* Fx(:,6) + ...    %   Fx23
                ((H_Vec(3)*H_Vec(1))/3 - H_Vec(1)*params.kz .*  kH_over_k2) .* Fx(:,7) + ...     %  Fx31
                ((H_Vec(3)*H_Vec(2))/3 - H_Vec(2)*params.kz .*  kH_over_k2) .* Fx(:,8) + ...     %  Fx32
                ((H_Vec(3)^2)/3 - H_Vec(3)*params.kz .* kH_over_k2) .* Fx(:,9) + ...                  %  Fx33
                Fx(:,10);                                                                                                                  %  Offset
        end   
    res = Res(:);   
    fprintf('+')
    end
else
    error('solving failed: undefined solving model');
end