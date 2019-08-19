function [dA_kt] = computeDistorsion(X_ft,Y_ftk,NMFparams)
%% Load Parameters
B = NMFparams.B;
k_max = NMFparams.k_max;
t_max = NMFparams.t_max;

%% Computing B-Distorsion
dA_kt = zeros(k_max,t_max);

for kk = 1:k_max
    if B==1
        dA_kt(kk,:) = sum( (X_ft .* log((X_ft./(squeeze(Y_ftk(:,:,kk))+eps))+eps)) - X_ft + squeeze(Y_ftk(:,:,kk)) ,1);
    elseif B==0
        dA_kt(kk,:) = sum( (X_ft./(squeeze(Y_ftk(:,:,kk))+eps)) - log((X_ft./(squeeze(Y_ftk(:,:,kk))+eps))+eps) - 1 ,1);
    else
        dA_kt(kk,:) = sum( (1/(B*(B-1))) * (X_ft.^B + (B-1)*(squeeze(Y_ftk(:,:,kk))+eps).^B- B*X_ft.*(squeeze(Y_ftk(:,:,kk))+eps).^(B-1)) ,1);
    end
end
end
