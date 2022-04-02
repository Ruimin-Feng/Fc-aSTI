function [ mae ] = compute_mae( chi_recon, chi_true,mask )


mae = sum(abs( chi_recon(mask==1) - chi_true(mask==1) )) / sum(mask(:));


end