function [AE_list,AE_mean,AE_std]=compute_AE(Vt,Vr,mask)
%========================
% This function calculates the angular error between Vt and Vr.
% Inputs:
% Vt: target PEV maps
% Vr: reconstructed PEV maps
%========================
       temp=acos(abs(dot(Vt,Vr,4)))/pi*180.*mask;
      AE_list=temp(mask==1);
      AE_mean=mean(AE_list);
      AE_std=std(AE_list);      
end