function [cfreq_amplitudes,miditobins,indval,minkmin]=computeCfreqFast(x,param)
% [cfreq_amplitudes,miditobins]=computeCfreqFast(x,param)
% Define the default parameters from the input signal
%
% Inputs:
% x  - signal
% param - input parameters
%
% Outputs:
% cfreq_amplitudes - amplitude "MIDI" spectrogram (mnotesXframes)
% miditobins - Relation MIDI note & MIDI bins (2Xmnotes)

% Inicializaciones
fs = param.fs;
longitud_frec = param.longitud_frec;
n_tramas = param.n_tramas;
midi_inc = param.midi_inc;
midi_min = param.midi_min;
midi_max = param.midi_max;
windowsize = param.windowsize;
salto = param.salto;

ventana = sqrt(hanning(windowsize,'periodic'));
% FFT Computation
X=sg(x,2*longitud_frec,fs,ventana,windowsize-salto);

miditobinskmin = zeros(1,(midi_max-midi_min+1)*midi_inc);
miditobinskmax = zeros(1,(midi_max-midi_min+1)*midi_inc);

for nota_midi=midi_min:midi_max
    
    for midi_interval = 0:midi_inc-1
        
        step_midi = 1/midi_inc;
        
        fmin=((2^(((nota_midi+midi_interval*step_midi-step_midi/2)-69)/12))*440);
        kmin=ceil(fmin/fs*(2*longitud_frec)+1);
        kmin=min(kmin,longitud_frec+1);
        
        fmax=((2^(((nota_midi+midi_interval*step_midi+step_midi/2)-69)/12))*440);
        kmax=fix(fmax/fs*(2*longitud_frec)+1);
        kmax=min(kmax,longitud_frec);
        
        miditobinskmin((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmin;
        miditobinskmax((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmax;
        
    end
end

if midi_inc==1
    miditobinskmax(miditobinskmax<miditobinskmin)=miditobinskmin(miditobinskmax<miditobinskmin);
else
    indval = (miditobinskmax>=miditobinskmin);
    miditobinskmin = miditobinskmin(indval);
    miditobinskmax = miditobinskmax(indval);
    minkmin = min(miditobinskmin);
    miditobinskmin = [1:minkmin-1 miditobinskmin];
    miditobinskmax = [1:minkmin-1 miditobinskmax];
end

miditobinskmax(end) = longitud_frec+1;
muestrasmidi = length(miditobinskmin);
miditobins=[miditobinskmin;miditobinskmax];

cfreq_amplitudes = zeros(n_tramas,muestrasmidi);
cfreq_amplitudes = cfreq_amplitudes';

for midi_index=1:muestrasmidi
    
    kmin=miditobinskmin(midi_index);
    kmax=miditobinskmax(midi_index);
    
    if (kmax-kmin)==0
        cfreq_amplitudes(midi_index,:) = abs(X(kmin,:));
    elseif (kmax-kmin)>0
        cfreq_amplitudes(midi_index,:) = sqrt(sum(abs(X(kmin:kmax,:)).^2,1)); % / (kmax-kmin+1));
    end
end

return;