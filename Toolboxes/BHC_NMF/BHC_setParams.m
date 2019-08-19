function [NMFparams,X_ft]=BHC_setParams(cfreq_amplitudes,param,NMFparams)
% NMF function parameters
%
% [NMFparams,X_ft]=BHC_setParams(cfreq_amplitudes,param,NMFparams)
%
% NMFparams.
%           NUM_MAX_ITER - Maximum number of NMF iterations (150)
%           ALPHA_A - Temporal smoothness value of A (50)
%           B - Beta divergence (0 itakura, 1 Kullback, 2 EUC) (1)
%           ETA_C - Update convergence speed for C (0.5)
%           ETA_A - Update convergence speed for A (0.5)
%           n_max - Excitation filter number (1)
%           m_max - number of partials (20)
% Input arguments: 
%	cfreq_amplitudes = spectral amplitude x time matrix
%   param = Sinusoidal model parameters
%   NMFparams* = Create if not exits-
%
% Output: 
%	NMFparams = NMF function parameters
%
% Julio Carabias / Francisco Rodriguez Diciembre 2011

% Si no existe NMFparams, lo creo
if nargin<3,
    NMFparams=struct();
end;

% Numero maximo de iteraciones
if ~isfield(NMFparams,'NUM_MAX_ITER')
    NMFparams.NUM_MAX_ITER=150;end;
% Valor de suavidad para los pesos
if ~isfield(NMFparams,'ALPHA_A')
    NMFparams.ALPHA_A=50;end;
% Valor sparseness penalty term
if ~isfield(NMFparams,'lambda')
    NMFparams.lambda=0;end;  
% Beta divergencia (0 itakura, 1 Kullback, 2 EUC)
if ~isfield(NMFparams,'B')
    NMFparams.B=1;end;
% Coeficente de retardo en la actualizacion de C
if ~isfield(NMFparams,'ETA_C')
    NMFparams.ETA_C=0.5;end;
% Coeficente de retardo en la actualizacion de A
if ~isfield(NMFparams,'ETA_A')
    NMFparams.ETA_A=0.5;end;  

% Numero de excitation filters
if ~isfield(NMFparams,'n_max')
NMFparams.n_max=1;end;

% Numero de parciales
if ~isfield(NMFparams,'m_max')
NMFparams.m_max=param.delta_max;
end;

% Numero de instrumentos
if ~isfield(NMFparams,'j_max')
NMFparams.j_max=1;end;

% Excitaciones con forma de ventana
if ~isfield(NMFparams,'wind_excit')    
    NMFparams.wind_excit = nmf_generateWindExcit(param);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of harmonics
NMFparams.harmonicos=unique(round(12*log2(1:NMFparams.m_max)));

X_ft = cfreq_amplitudes(:,param.midi_min:param.midi_max)';

NMFparams.f_max=size(X_ft,1);
NMFparams.t_max=size(X_ft,2); 
NMFparams.p_max=param.midi_max-param.midi_min+1;
NMFparams.p_min = param.midi_min;
