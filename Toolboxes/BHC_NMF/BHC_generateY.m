function [Y_ft]=BHC_generateY(NMFparams,S_fpj,A_ptj)
% Create the estimation Y
%
% [Y_ft]=nmf_generateY(NMFparams,S_fpj,A_ptj)
%
% Y_ft = sum_p sum_j ( S_fpj A_ptj )
%
% Julio Carabias / Francisco Rodriguez -- March 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_max = NMFparams.f_max;
t_max = NMFparams.t_max;
j_max = NMFparams.j_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_ft=zeros(f_max,t_max);
for jj=1:j_max,
    Y_ft = Y_ft + (S_fpj(:,:,jj) * A_ptj(:,:,jj));
end;

return;
