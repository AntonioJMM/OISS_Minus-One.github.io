function [C_pmj_INC] = BHC_INCmodel(C_pmj,incremento)
% function [C_pmj_INC] = BHC_INCmodel(C_pmj,incremento)
% 
% Generate an extended BHC instrument base models increasing the frequency
% resolution
% 
% Inputs:
% - C_pmj: Original model (midi_inc = 1)
% - incremento: N times more freq resolution (midi_inc)
%
C_pmj_INC = zeros(size(C_pmj,1)*incremento,size(C_pmj,2),size(C_pmj,3));

for q=1:size(C_pmj,1)
    if q==1
        C_pmj_INC(1:floor(incremento/2),:) = repmat(C_pmj(q,:),floor(incremento/2),1);
    else
        C_pmj_INC((q-1)*incremento-(incremento/2-1):((q-1)*incremento+(incremento/2)),:) = repmat(C_pmj(q,:),incremento,1);     
    end
end

return;