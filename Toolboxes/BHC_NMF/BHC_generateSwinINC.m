function [S_fpj]=BHC_generateSwinINC(NMFparams,C_pmj)
% Create the basis S
%
% [S_fpj]=BHC_generateSwinINC(NMFparams,C_pmj)
%
% S_fpj = sum_n sum_m ( H_jf G_nmj K_npj )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wind_excit = NMFparams.wind_excit;
f_max = 401;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;

if f_max>114
    harmonicos_pitches = NMFparams.harmonicos_pitches;
else
    harmonicos = NMFparams.harmonicos;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_fpj=zeros(f_max,p_max,j_max);

for p=1:p_max
    if f_max>114
        deltas=harmonicos_pitches(p,:);
        deltas = deltas(deltas~=0);
        deltas=deltas(deltas<=f_max);
        num_deltas=length(deltas);
    else
        deltas=harmonicos+p;
        deltas = deltas(deltas~=0);
        deltas=deltas(deltas<=p_max);
        num_deltas=length(deltas);
    end
    
    for jj=1:j_max
        S_fpj(:,p,jj) = C_pmj(p,1:num_deltas,jj) * wind_excit(deltas,:);
    end
end

return;
