function [temp,AE_mean,AE_std]=compute_AE(Vt,Vr,mask)
     norm_Vt=sqrt(Vt(:,:,:,1).^2+Vt(:,:,:,2).^2+Vt(:,:,:,3).^2);
     norm_Vr=sqrt(Vr(:,:,:,1).^2+Vr(:,:,:,2).^2+Vr(:,:,:,3).^2);
      temp=real(acos(abs(dot(Vt,Vr,4)./(norm_Vt.*norm_Vr))))/pi*180.*mask;
      AE_list=temp(mask==1);
      AE_mean=mean(AE_list);
      AE_std=std(AE_list);      
end