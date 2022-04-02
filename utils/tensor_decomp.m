function [symm,asymm]=tensor_decomp(sti)
SS=size(sti);
sti=reshape(sti,[SS(1),SS(2),SS(3),3,3]);
sti_t=permute(sti,[1,2,3,5,4]);
symm=reshape(0.5*(sti+sti_t),[SS(1),SS(2),SS(3),9]);
asymm=reshape(0.5*(sti-sti_t),[SS(1),SS(2),SS(3),9]);

end