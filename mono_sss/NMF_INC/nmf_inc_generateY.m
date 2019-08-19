function [Y_ft]=nmf_inc_generateY(NMFparams,S_fpj,A_ptj,eS_pj,eA_t)

f_max = NMFparams.f_max;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;

S_fpj = S_fpj.*permute(repmat(eS_pj,1,1,f_max),[3 1 2]);
A_ptj = A_ptj.*repmat(eA_t,p_max,1,j_max);
Y_ft = tprod(S_fpj,[1 -1 -2],A_ptj,[-1 2 -2]);

return;
