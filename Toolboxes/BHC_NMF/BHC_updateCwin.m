function [C_pmj]=BHC_updateCwin(NMFparams,C_pmj,A_ptj,X_ft,Y_ft)
% Update basis S with
%
% [C_pmj]=BHC_updateCwin(NMFparams,C_pmj,A_ptj,X_ft,Y_ft)
%
%                   sum_t sum_f ( A_ptj W_pf(deltas,:) Y_ft^(B-2) X_ft )
% C_pmj <- C_pmj * -----------------------------------------------------------
%                   sum_t sum_f ( A_ptj W_pf(deltas,:) Y_ft^(B-1) )
%
% Finally apply a filter to smooth the instrument envelope
%
% Julio Carabias / Francisco Rodriguez -- March 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = NMFparams.B;
ETA_C = NMFparams.ETA_C;


harmonicos = NMFparams.harmonicos;
f_max = NMFparams.f_max;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;

wind_excit = NMFparams.wind_excit;
wind_excit_trans = wind_excit';
sparse_pen = NMFparams.lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Valores de actualizacion
CactN=C_pmj.*0;
CactD=C_pmj.*0;

%% Computo actualizacion segun gradiente
valor1=X_ft .* (Y_ft.^(B-2));
valor1(~isfinite(valor1))=0;
valor2=Y_ft.^(B-1);
valor2(~isfinite(valor2))=0;

for p=1:p_max,
    % Frecuencias de cada parcial del pitch actual
    deltas=intersect([1:f_max],harmonicos+p);
    num_deltas=length(deltas);
    % Para cada parcial
    for m=1:num_deltas,
        % f es la frecuencia del parcial actual
        f=deltas(m);
        for jj=1:j_max
            A_wind = wind_excit_trans(:,f) * A_ptj(p,:,jj);
            
            CactN(p,m,jj)=sum(A_wind(f,:) .* valor1(f,:));
            CactD(p,m,jj)=sum(A_wind(f,:) .* valor2(f,:));
        end; % Fin j
    end; % Fin m
end; % Fin p

C_pmj = C_pmj .* (CactN./(CactD+sparse_pen)).^ETA_C;
C_pmj(~isfinite(C_pmj))=0;

return;
