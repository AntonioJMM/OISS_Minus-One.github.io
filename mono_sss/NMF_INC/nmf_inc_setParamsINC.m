function NMFparams = nmf_inc_setParamsINC(param,NMFparams)
% NMF function parameters
%
% NMFparams.
%           NUM_MAX_ITER - Maximum number of NMF iterations (150)
%           ALPHA_A - Temporal smoothness value of A (50)
%           B - Beta divergence (0 itakura, 1 Kullback, 2 EUC) (1)
%           ETA_C - Update convergence speed for C (0.5)
%           ETA_A - Update convergence speed for A (0.5)
%           ETA_M - Update convergence speed for M (0.5)
%           n_max - Excitation filter number (1)
%           m_max - number of partials (20)
%           j_max - number of instruments
%           wind_excit - window transform
%           NMFparams.harmonicos_pitches
%           NMFparams.harmonicos
%           NMFparams.fftparam=param;
% Input arguments:
%	cfreq_amplitudes = spectral amplitude x time matrix
%   param = Sinusoidal model parameters
%   NMFparams* = Create if not exits-
%
% Output:
%	NMFparams = NMF function parameters

% Numero maximo de iteraciones
if ~isfield(NMFparams,'NUM_MAX_ITER')
    NMFparams.NUM_MAX_ITER=150;
end
% Valor de suavidad para los pesos
if ~isfield(NMFparams,'ALPHA_A')
    NMFparams.ALPHA_A=50;
end
% Valor sparseness penalty term
if ~isfield(NMFparams,'lambda')
    NMFparams.lambda=0;
end
% Beta divergencia (0 itakura, 1 Kullback, 2 EUC)
if ~isfield(NMFparams,'B')
    NMFparams.B=1;
end
% Coeficente de retardo en la actualizacion de C
if ~isfield(NMFparams,'ETA_C')
    NMFparams.ETA_C=0.5;
end
% Coeficente de retardo en la actualizacion de A
if ~isfield(NMFparams,'ETA_A')
    NMFparams.ETA_A=0.5;
end
% Coeficente de retardo en la actualizacion de M
if ~isfield(NMFparams,'ETA_M')
    NMFparams.ETA_M=0.5;
end
% Numero de excitation filters
if ~isfield(NMFparams,'n_max')
    NMFparams.n_max=1;
end
% Numero de parciales
if ~isfield(NMFparams,'m_max')
    NMFparams.m_max=param.delta_max;
end
% Numero de instrumentos
if ~isfield(NMFparams,'j_max')
    NMFparams.j_max=1;
end
% Excitaciones con forma de ventana
if ~isfield(NMFparams,'wind_excit')
    NMFparams.wind_excit=nmf_generateWindExcitINC(param);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NMFparams.wind_width = getWindWidth(NMFparams.wind_excit);

NMFparams.p_max=(param.midi_max-param.midi_min+1)*param.midi_inc;
NMFparams.p_min = 1;

% Number of harmonics
%NMFparams.harmonicos=unique(round(12*log2(1:NMFparams.m_max)));
NMFparams.flog_first = param.flog_first;
harmonicos = zeros(NMFparams.p_max,param.delta_max);
harmonicos_pitches = zeros(NMFparams.p_max,NMFparams.m_max);
NMFparams.pitches=param.midi_min:1/param.midi_inc:param.midi_max+1-1/param.midi_inc;

for p=1:NMFparams.p_max
    for d=1:NMFparams.m_max
        midi_INC_of_p=p;
        midi_act_of_p = param.midi_min+((midi_INC_of_p-1)/param.midi_inc);
        frec_act_of_p = ((2^((midi_act_of_p-69)/12))*440);
        
        frec_act = d*frec_act_of_p;
        midi_act = 12*log2(frec_act/440)+69;
        if ((round(midi_act*param.midi_inc)/param.midi_inc)<=(NMFparams.pitches(end)))
            harmonicos_pitches(p,d) = (((round(midi_act*param.midi_inc)/param.midi_inc)-param.midi_min)*param.midi_inc)+1;
        end
        bin_act = (frec_act*param.longitud_frec*2/param.fs);
        if ((round(bin_act)+1)<=(param.longitud_frec+1))
            harmonicos(p,d) = param.bins2flog(round(bin_act)+1);
        end
    end
end

NMFparams.harmonicos_pitches = harmonicos_pitches;
NMFparams.harmonicos=unique(round(12*log2(1:NMFparams.m_max)));
NMFparams.fftparams=param;
NMFparams.f_max = 401;

return;

%--------------------------------------------------------------

function [win_width] = getWindWidth(wind_excit)

aux = 20*log10(wind_excit);
aux=aux>-30;
win_width = sum(aux,2);

return;