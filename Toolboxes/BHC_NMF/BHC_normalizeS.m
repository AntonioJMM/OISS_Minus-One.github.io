function [S_fpj, C_pmj, A_ptj]=BHC_normalizeS(NMFparams,S_fpj,C_pmj,A_ptj)
%% [S_fpj, C_pmj, A_ptj]=BHC_normalizeS(NMFparams,S_fpj,C_pmj,A_ptj)
% Normaliza las bases para norm1 = 1.
% Julio Carabias / Francisco Rodriguez -- December 2011

p_max = NMFparams.p_max;
j_max = NMFparams.j_max;

% Normalizo S sobre C y A
for jj=1:j_max,
    for p=1:p_max,
        escala=sqrt(sum(S_fpj(:,p,jj).^2));
        
        % Escala S
        S_fpj(:,p,jj)=S_fpj(:,p,jj)./escala;
        S_fpj(~isfinite(S_fpj))=0;

        % Escala C
        C_pmj(p,:,jj)=C_pmj(p,:,jj)./escala;
        C_pmj(~isfinite(C_pmj))=0;
        
        % Traslada factor de escala a A
        A_ptj(p,:,jj)=A_ptj(p,:,jj).*escala;
    end;
end;

return;
