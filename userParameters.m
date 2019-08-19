% SOURCE SEPARATION - USER PARAMETERS

%% USER --> SET RELATIVE PATHS
% Paths from BD URMP
param.ruta_bbdd    = [pwd '/BBDD/URMP/'];
param.ruta_modelos = [pwd '/BBDD/INSTRUMENT_MODELS/'];
param.outpath      = [pwd '/OUTPUTS/'];
if exist(param.outpath,'dir')==0, mkdir(param.outpath); end

%% USER --> FFT INPUT PARAMETERS
param.midi_inc=4;
param.midi_min=24;
param.midi_max=137;
param.t_trama=0.128;
param.t_salto=0.032;
param.longitud_frec=8192;

%% USER --> NMF INPUT PARAMETERS
NMFparams.ETA_A = 0.5;
NMFparams.ETA_M = 0.1;
NMFparams.NUM_MAX_ITER = 50;
NMFparams.ALPHA_A = 0;      % Temporal smoothness for the gains
NMFparams.B = 1.3;          % Beta divergencia (0 itakura, 1 Kullback, 2 EUC)
NMFparams.m_max = 20;       % Number of harmonics
NMFparams.n_max = 0;        % Number of excitations
