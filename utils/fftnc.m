function im = fftnc(d)
% im = fftnc(d)
%
% fftnc perfomrs a centered fftn
%
im = fftshift(fftn(fftshift(d)));
end

