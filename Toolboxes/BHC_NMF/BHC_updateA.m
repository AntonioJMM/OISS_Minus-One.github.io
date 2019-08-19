function [A_ptj]=BHC_updateA(NMFparams,A_ptj,S_fpj,X_ft,Y_ft)
% Update gains A with
%
% [A_ptj]=nmf_updateA(NMFparams,A_ptj,S_fpj,X_ft,Y_ft)
%
%                     sum_f ( S_fpj Y_ft^(B-2) X_ft )
% A_ptj <- A_ptj * --------------------------------------
%                       sum_f ( S_fpj Y_ft^(B-1) )
%
% if B==1 (divergence), applies Virtanen temporal statitiscal constrains
%
%                    2a/A_ptj * sum_f ( S_fpj Y_ft^(B-2) X_ft )
% A_ptj <- A_ptj * ----------------------------------------------------
%                    a(Z_ptj + Z_p(t-1)j) * sum_f ( S_fpj Y_ft^(B-1) )
%
% where
%       |   1 / A_p1j             : t==1
% Z  <- |   2/(A_ptj + A_p(t-1)j) : 1 < t < T+1
%       |   1/A_ptj               : t== T+1
%
% Julio Carabias / Francisco Rodriguez -- March 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B       = NMFparams.B;
ETA_A   = NMFparams.ETA_A;
ALPHA_A = NMFparams.ALPHA_A;

t_max = NMFparams.t_max;
p_max = NMFparams.p_max;
j_max = NMFparams.j_max;

sparse_pen = NMFparams.lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Temporal Statistical Smoothness (Virtanen)
ApenN=A_ptj.*0;
ApenD=A_ptj.*0;
if B==1 && ALPHA_A>0,
    temp_a=ALPHA_A;
    temp_b=0;
    
    % Genero factor de penalizacion
    Z_ptj=zeros(p_max,t_max+1,j_max);
    
    for jj=1:j_max,
        for p=1:p_max,
            Z_ptj(p,1,jj)=1./(A_ptj(p,1,jj)+temp_b);
            Z_ptj(~isfinite(Z_ptj))=0;
            for t=2:t_max,
                Z_ptj(p,t,jj)=2./(A_ptj(p,t,jj) + A_ptj(p,t-1,jj));
                Z_ptj(~isfinite(Z_ptj))=0;
            end;
            Z_ptj(p,t_max+1,jj)=1./(A_ptj(p,t_max,jj));
            Z_ptj(~isfinite(Z_ptj))=0;
        end; % Fin p
    end; % Fin j
    
    % Valores de actualizacion
    ApenN=(2*temp_a) ./ (A_ptj);
    ApenN(~isfinite(ApenN))=0;
    ApenD=temp_a .* (Z_ptj(:,1:end-1,:) + Z_ptj(:,2:end,:));
end;

%% Factor Sparsity
if NMFparams.lambda>0,
    sparse_pen = zeros(p_max,t_max,j_max);
    for jj=1:j_max,
        sparse_factor = NMFparams.lambda * mean(sum(A_ptj(:,:,jj),2));
        sparse_pen(:,:,jj)=sparse_factor;
    end;
end;
        
%% Classical NMF gains update

% Valores de actualizacion
AactN=A_ptj.*0;
AactD=A_ptj.*0;

% Valores de X*Y^B-2 y Y^B-1 
valor1=X_ft .* (Y_ft.^(B-2));
valor1(~isfinite(valor1))=0;
valor2=Y_ft.^(B-1);
valor2(~isfinite(valor2))=0;

% Resto del sumatorio
for p=1:p_max,
    for jj=1:j_max        
        AactN(p,:,jj)=sum(bsxfun(@times,valor1,S_fpj(:,p,jj)),1);
        AactD(p,:,jj)=sum(bsxfun(@times,valor2,S_fpj(:,p,jj)),1);
    end; % Fin j
end; % Fin p

% Resultados
AactN=AactN + ApenN;
AactD=AactD + ApenD + sparse_pen;
A_ptj = A_ptj .* (AactN./(AactD)).^ETA_A;
A_ptj(~isfinite(A_ptj))=0;

return;