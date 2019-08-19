function [param]=getParametrosMusica(param,nsamples,fs)
% [param]=getParametrosMusica(param,nsamples,fs)
% Obtain the parameters associated to the time frequency representation
%
% Input arguments:
%	param --> see userparameter.m to set initial values
%   nsamples -->number of samples of the input signal
%   fs --> sampling frequnecy
%
% Output:
%   param --> Initialized parameters fo the TF representation

param.fs=fs;
param.ts=1/param.fs;

% Parameter Definition
if ~isfield(param,'midi_min')
    param.midi_min=24;
end
if ~isfield(param,'midi_inc')
    param.midi_inc=0;
end
if ~isfield(param,'delta_max')
    param.delta_max=9;
end
if ~isfield(param,'npolos')
    param.npolos=1;
end
if ~isfield(param,'t_trama')
    param.t_trama=0.128;
end
if ~isfield(param,'t_salto')
    param.t_salto=0.032;
end
if ~isfield(param,'longitud_frec')
    param.longitud_frec=8192;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.windowsize=round(param.t_trama*fs/2)*2;
param.salto=round(param.t_salto*fs/2)*2;

% Actualizamos los valores de tiempo
param.t_trama=param.windowsize/fs;
param.n_tramas=fix((nsamples-param.windowsize)/param.salto)+1;

% Configuracion del modelado sinusoidal
param.midi_max=ceil(12*log2((fs/2)/440)+69);
param.umbral_perceptual=.05;
param.numtonos_maximo=100;

%Cálculo de indices válidos
miditobinskmin = zeros(1,(param.midi_max-param.midi_min+1)*param.midi_inc);
miditobinskmax = zeros(1,(param.midi_max-param.midi_min+1)*param.midi_inc);
for nota_midi=param.midi_min:param.midi_max
    for midi_interval = 0:param.midi_inc-1
        
        step_midi = 1/param.midi_inc;
        
        fmin=((2^(((nota_midi+midi_interval*step_midi-step_midi/2)-69)/12))*440);
        kmin=ceil(fmin/param.fs*(2*param.longitud_frec)+1);
        kmin=min(kmin,param.longitud_frec+1);
        
        fmax=((2^(((nota_midi+midi_interval*step_midi+step_midi/2)-69)/12))*440);
        kmax=fix(fmax/param.fs*(2*param.longitud_frec)+1);
        kmax=min(kmax,param.longitud_frec);
        
        miditobinskmin((nota_midi-param.midi_min)*param.midi_inc+midi_interval+1) = kmin;
        miditobinskmax((nota_midi-param.midi_min)*param.midi_inc+midi_interval+1) = kmax;
        
    end
end

indval = (miditobinskmax>=miditobinskmin);
miditobinskmin = miditobinskmin(indval);
miditobinskmax = miditobinskmax(indval);
minkmin = min(miditobinskmin);
miditobinskmin = [1:minkmin-1 miditobinskmin];
miditobinskmax = [1:minkmin-1 miditobinskmax];
miditobinskmax(end) = param.longitud_frec+1;
muestrasmidi_inc = length(miditobinskmin);
midi_inc2bins=[miditobinskmin;miditobinskmax];

midi_cont=1;
param.bins2flog = zeros(1,param.longitud_frec+1);
for bin=1:param.longitud_frec+1
    if(midi_inc2bins(2,midi_cont)<bin)
        midi_cont = midi_cont+1;
    end
    param.bins2flog(bin) = midi_cont;
end

param.flog2bins = midi_inc2bins;
param.indval = indval;
param.muestrasflog = muestrasmidi_inc;
param.flog_first = minkmin;

return