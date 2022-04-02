function out=add_noise(in_data,SNR,mask)
     NOISE=randn(size(in_data));
     NOISE=(NOISE-mean(NOISE(:)))/std(NOISE(:));
     signal_power=norm(in_data(mask==1))^2/sum(mask(:));
     noise_var=signal_power/(10^(SNR/10));
     NOISE=sqrt(noise_var)*NOISE;
     out=(in_data+NOISE).*mask;    

end