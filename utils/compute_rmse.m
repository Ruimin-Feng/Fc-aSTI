function [ rmse ] = compute_rmse( chi_recon, chi_true )
%=========================
% This function calculates RMSE
% Inputs:
% chi_recon: reconstructed maps
% chi_true: ground truth maps
%=========================

rmse = 100 * norm( chi_recon(:) - chi_true(:) ) / norm(chi_true(:));


end

