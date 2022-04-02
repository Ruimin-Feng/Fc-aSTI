function im = ifftnc(d)
% im = ifftnc(d)
%
% ifftnc perfomrs a centered ifftn
%
im = fftshift(ifftn(fftshift(d)));
end
