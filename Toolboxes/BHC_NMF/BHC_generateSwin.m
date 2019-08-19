function [S_fpj]=BHC_generateSwin(NMFparams,C_pmj)
% Create the basis S
%
% [S_fpj]=nmf_generateS(NMFparams,H_jf,G_nmj,K_npj)
%
% S_fpj = sum_n sum_m ( H_jf G_nmj K_npj )
%
% Julio Carabias / Francisco Rodriguez. March 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
harmonicos = NMFparams.harmonicos;
wind_excit = NMFparams.wind_excit;

f_max = NMFparams.f_max;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_fpj=zeros(f_max,p_max,j_max);

for p=1:p_max,
    
    deltas = harmonicos+p;
    deltas = deltas(deltas~=0);
    deltas = deltas(deltas<=p_max);
    num_deltas = length(deltas);
    
    for jj=1:j_max,
        S_fpj(:,p,jj) = C_pmj(p,1:num_deltas,jj) * wind_excit(deltas,:);
    end;
    
end;

return;
