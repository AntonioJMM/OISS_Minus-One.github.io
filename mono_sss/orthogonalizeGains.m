function [A_1t] = orthogonalizeGains(X_ft,W_f,NMFparams)
% [A_1t] = orthogonalizeGains(X_ft,S_fpj,E_1pj)
% Genera la proyeccion ortogonal para una combinacion concreta en todos los
% tiempo de ocurrencia.
%
% Inputs:
% NMFparams - NMF default params 
% X_ft      - Signal
% W_f     - Bases
%
% Outputs:
% A_1t      - Projected gains
%
%% Load Parameters
t_max = size(X_ft,2);
B     = NMFparams.B;

%% Generar parámetros de salida
A_1t = zeros(1,t_max);

%% Estimacion de ganancias
A_1t = sum( X_ft.*repmat(W_f.^(B-1),1,t_max) ,1) ...
    / sum( W_f.^B + eps ,1);

end
