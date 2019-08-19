function [E_kpj] = nmf_inc_updateE(X_ft,NMFparams,S_fpj,E_kpj,A_kt)

%% Load Parameters
NUM_MAX_ITER = 50;
B = NMFparams.B;

f_max = NMFparams.f_max;
t_max = NMFparams.t_max;
p_max = NMFparams.p_max; % pitch
j_max = NMFparams.j_max; % instrumentos
k_max = NMFparams.k_max; % combinaciones
s_max = 1;

%% Normalizamos la senal de entrada
if (B==1)
    norma_ft = (sum(sum(X_ft.^B))^(1/B));
    Xn_ft = X_ft/norma_ft;
elseif (B==0)
    norma_ft = size(X_ft,1)*size(X_ft,2);
    Xn_ft = X_ft/norma_ft;
else
    norma_ft = (sum(sum(X_ft.^B))^(1/B)) / ((B*(B-1)).^(1/B));
    Xn_ft = X_ft/norma_ft;
end

%% Inicializo M
M_js = ones(j_max,s_max);

%% Normalizaciones
eS_pj = sqrt(sum(S_fpj.^2,1));
Sn_fpj = S_fpj./repmat(eS_pj,f_max,1,1);
Sn_fpj(~isfinite(Sn_fpj))=0;

eS_pj = ones(p_max,j_max);

eE_k = sqrt(sum(sum(E_kpj.^2,2),3));
En_kpj = E_kpj./repmat(eE_k,1,p_max,j_max);
En_kpj(~isfinite(En_kpj))=0;

An_kt = bsxfun(@times, A_kt, eE_k); 

eA_t = sqrt(sum(An_kt.^2,1));
An_kt = An_kt./repmat(eA_t,k_max,1);
An_kt(~isfinite(An_kt))=0;

eA_t = ones(1,t_max);

An_ptj = tprod(En_kpj,[-1 1 3],An_kt,[-1 2]);

[Y_ft]=nmf_inc_generateY(NMFparams,Sn_fpj,An_ptj,eS_pj,eA_t);
if (B==1)
    norma_Yft = (sum(sum(Y_ft.^B,2))).^(1/B);
elseif (B==0)
    norma_Yft = size(Y_ft,1)*size(Y_ft,2);
else
    norma_Yft = (sum(sum(Y_ft.^B,2)).^(1/B)) ./ ((B*(B-1))^(1/B));
end

eS_pj = eS_pj/norma_Yft;
[Y_ft]=nmf_inc_generateY(NMFparams,Sn_fpj,An_ptj,eS_pj,eA_t);


%% Costes
D_B = zeros(NUM_MAX_ITER+1,1);

if B==1
    D_B(1) = sum(sum((Xn_ft .* log((Xn_ft./(Y_ft+eps))+eps)) - Xn_ft + Y_ft));
elseif B==0
    D_B(1) = sum(sum((Xn_ft./(Y_ft+eps)) - log((Xn_ft./(Y_ft+eps))+eps) - 1));
else
    D_B(1) = sum(sum( (1/(B*(B-1))) * (Xn_ft.^B + (B-1)*(Y_ft+eps).^B- B*Xn_ft.*(Y_ft+eps).^(B-1)) ));
end

%% Iteraciones
for ite = 1:NUM_MAX_ITER
    %% States
    Y_fts_cb2 = 0*Y_ft;
    Y_fts_cb2((Y_ft~=0)) = Y_ft((Y_ft~=0)).^(B-2);
    Y_fts_cb1 = 0*Y_ft;
    Y_fts_cb1((Y_ft~=0)) = Y_ft((Y_ft~=0)).^(B-1);
    
    Sne_fpj = Sn_fpj.*permute(repmat(eS_pj,1,1,f_max),[3 1 2]);
    
    KL_num = tensor_mult(Sne_fpj,(Y_fts_cb2 .* Xn_ft),1,1);
    KL_num = permute(tensor_mult(KL_num,(An_kt * diag(eA_t)),3,2),[3 1 2]);
    
    KL_den = tensor_mult(Sne_fpj,Y_fts_cb1,1,1);
    KL_den = permute(tensor_mult(KL_den,(An_kt * diag(eA_t)),3,2),[3 1 2]);
    
    %% Actualization
    delta = KL_num ./ KL_den;
    delta(~isfinite(delta)) = 0;
    
    En_kpj = En_kpj .* delta;
    
    % Normalizacion
    eE_k = sqrt(sum(sum(En_kpj.^2,2),3));
    En_kpj = En_kpj./repmat(eE_k,1,p_max,j_max);
    En_kpj(~isfinite(En_kpj))=0;
    
    An_kt = bsxfun(@times, An_kt, eE_k);
    
    eA_t_temp = sqrt(sum(An_kt.^2,1));
    eA_t_temp(~isfinite(eA_t_temp)) = 0;
    rep=repmat(eA_t_temp,k_max,1);
    
    An_kt(rep~=0) = An_kt(rep~=0)./rep(rep~=0);
    An_kt(~isfinite(An_kt))=0;
    
    eA_t = eA_t .* eA_t_temp;
    
    An_ptj = tprod(En_kpj,[-1 1 3],An_kt,[-1 2]);
    
    %% Reconstruction y normalizacion
    [Y_ft]=nmf_inc_generateY(NMFparams,Sn_fpj,An_ptj,eS_pj,eA_t);
    if (B==1)
        norma_Yft = (sum(sum(Y_ft.^B,2))).^(1/B);
    elseif (B==0)
        norma_Yft = size(Y_ft,1)*size(Y_ft,2);
    else
        norma_Yft = (sum(sum(Y_ft.^B,2)).^(1/B)) ./ ((B*(B-1))^(1/B));
    end
    
    eS_pj = eS_pj/norma_Yft;
    [Y_ft]=nmf_inc_generateY(NMFparams,Sn_fpj,An_ptj,eS_pj,eA_t);
    
    %% Costs
    if B==1
        D_B(ite+1) = sum(sum((Xn_ft .* log((Xn_ft./(Y_ft+eps))+eps)) - Xn_ft + Y_ft));
    elseif B==0
        D_B(ite+1) = sum(sum((Xn_ft./(Y_ft+eps)) - log((Xn_ft./(Y_ft+eps))+eps) - 1));
    else
        D_B(ite+1) = sum(sum( (1/(B*(B-1))) * (Xn_ft.^B + (B-1)*(Y_ft+eps).^B- B*Xn_ft.*(Y_ft+eps).^(B-1)) ));
    end
    
end

E_kpj = En_kpj;
return;