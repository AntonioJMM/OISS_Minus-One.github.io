function [NMF_inic,NMFparams,param] = alignmentStage(Xr_ft,NMF_inic,NMFparams,param)
% [NMF_inic,NMFparams,param] = alignmentStage(Xr_ft,NMF_inic,NMFparams,param)
% Compute the alignment stage
%
% Input arguments:
%   Xr_ft --> Input signal
%   NMF_inic --> Initialized NMF parameters
%   NMFparams --> store the NMF parameters
%	param --> see userparameter.m to set initial values
% 
% Output:
%   NMF_inic
%   NMFparams
%	param

%% Performing projection of gains
fprintf('\n\tPerforming projection of gains ... ');
B = NMFparams.B;
f_max = NMFparams.f_max;
E_kpj = NMF_inic.E_kpj;
S_fpj = NMF_inic.S_fpj;

if (B==1)
    norma_ft = (sum(sum(Xr_ft.^B))^(1/B));
    Xn_ft = Xr_ft/norma_ft;
elseif (B==0)
    norma_ft = size(Xr_ft,1)*size(Xr_ft,2);
    Xn_ft = Xr_ft/norma_ft;
else
    norma_ft = (sum(sum(Xr_ft.^B))^(1/B)) / ((B*(B-1)).^(1/B));
    Xn_ft = Xr_ft/norma_ft;
end

eS_pj = sqrt(sum(S_fpj.^2,1));
Sn_fpj = S_fpj./repmat(eS_pj,f_max,1,1);
Sn_fpj(~isfinite(Sn_fpj))=0;

W_fk = tprod(Sn_fpj,[1 -1 -2],E_kpj,[2 -1 -2]);

if (B==1)
    eW_k = (sum(W_fk.^B,1)).^(1/B);
elseif (B==0)
    eW_k = size(W_fk,1);
else
    eW_k = ((sum(W_fk.^B,1)).^(1/B)) ./ ((B*(B-1)).^(1/B));
end

Wn_fk = W_fk./repmat(eW_k,f_max,1);
Wn_fk(~isfinite(Wn_fk))=0;

A_kt = zeros(NMFparams.k_max,size(Xr_ft,2));
for kk = 2:NMFparams.k_max
    [A_kt(kk,:)] = orthogonalizeGains(Xn_ft,Wn_fk(:,kk),NMFparams);
end
fprintf('Done');

%% Aligment Process
fprintf('\n\tAligning gains ... ');
NMFparams.t_max = size(Xr_ft,2);

Y_ftk = tprod(Wn_fk,[1 3],A_kt,[3 2]);
[dA_kt] = computeDistorsion(Xn_ft,Y_ftk,NMFparams);
clear Y_ftk;

[NMF_inic.A_kt] = performDTW_function(dA_kt,param);
fprintf('Done');

return