function [cfreq_amplitudes,miditobins]=computeCfreqSeq(x,param,barra)
% [cfreq_amplitudes,miditobins]=computeCfreqSeq(x,param,barra)
% Define the default parameters from the input signal
%
% Inputs:
% x  - signal
% param - input parameters 
% barra - show process bar
%
% Outputs:
% cfreq_amplitudes - amplitude "MIDI" spectrogram (mnotesXframes)
% miditobins - Relation MIDI note & MIDI bins (2Xmnotes)
%
% Julio Carabias y Francisco Rodriguez. Fall 2012

% Inicializaciones
fs = param.fs;
longitud_frec = param.longitud_frec;
n_tramas = param.n_tramas;
midi_inc = param.midi_inc;
midi_min = param.midi_min;
midi_max = param.midi_max;

% Genero escala MIDI
miditobinskmin = zeros(1,(midi_max-midi_min+1)*midi_inc);
miditobinskmax = zeros(1,(midi_max-midi_min+1)*midi_inc);

for nota_midi=midi_min:midi_max,

    for midi_interval = 0:midi_inc-1,
        
        step_midi = 1/midi_inc;

        fmin=((2^(((nota_midi+midi_interval*step_midi-step_midi/2)-69)/12))*440);
        kmin=ceil(fmin/fs*(2*longitud_frec)+1);
        kmin=min(kmin,longitud_frec+1);

        fmax=((2^(((nota_midi+midi_interval*step_midi+step_midi/2)-69)/12))*440);
        kmax=fix(fmax/fs*(2*longitud_frec)+1);
        kmax=min(kmax,longitud_frec);

        miditobinskmin((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmin;
        miditobinskmax((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmax;
        
    end;

end;

if midi_inc==1,
    miditobinskmax(miditobinskmax<miditobinskmin)=miditobinskmin(miditobinskmax<miditobinskmin);
else   
    indval = (miditobinskmax>=miditobinskmin);
    miditobinskmin = miditobinskmin(indval);
    miditobinskmax = miditobinskmax(indval);
    minkmin = min(miditobinskmin); % Lineal hasta el primer midi
    miditobinskmin = [1:minkmin-1 miditobinskmin];
    miditobinskmax = [1:minkmin-1 miditobinskmax];
end;

miditobinskmax(end) = longitud_frec+1;
muestrasmidi = length(miditobinskmin);
miditobins=[miditobinskmin;miditobinskmax];

% Inicializacion de variables de salida
cfreq_amplitudes=zeros(n_tramas,muestrasmidi);
cfreq_amplitudes = cfreq_amplitudes';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculo el spectrograma secuencialmente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowsize = param.windowsize;
salto = param.salto;
nfft=2*longitud_frec;
win = sqrt(hanning(windowsize,'periodic'));  % Quitado SQRT a la ventana. FCO. 24/11/11
noverlap = windowsize-salto;

if ~any(any(imag(x)))    % x purely real
    if rem(nfft,2),    % nfft odd
        select = [1:(nfft+1)/2];
    else
        select = [1:nfft/2+1];
    end
else
    select = 1:nfft;
end

nx = length(x);
nwind = length(win);
if nx < nwind    % zero-pad x if it has length less than the win length
    x(end+1:nwind)=0;  nx=nwind;
end
x = x(:);     % make a column vector for ease later
win = win(:); % be consistent with data set

zp = zeros(nfft-nwind,1);
nwind2 = (nwind- mod(nwind,2))/2;
wzp = [win(nwind2+1:nwind);zp;win(1:nwind2)]; % !____!

%nframes = n_tramas;
nframes = fix((nx-noverlap)/(nwind-noverlap));
frame_index = 1 + (0:(nframes-1))*(nwind-noverlap);
if length(x)<(nwind+frame_index(nframes)-1)
    x(end+1:nwind+frame_index(nframes)-1) = 0;   % zero-pad x
end

% La puta barra
if barra,
    h = waitbar(0,'Please wait...','Name','Computing Cfreq ...');
end;

% X = zeros(length(select),nframes);
for t=1:nframes,
    % Actualizo la barra
    if barra && mod(t,10)==0,
        waitbar(t/n_tramas,h,['Trama ' num2str(t) ' de ' num2str(n_tramas)]);
    end;
    
    xframe = x(frame_index(t):frame_index(t)+nwind-1); % extract frame of input data
    xzp = [xframe(nwind2+1:nwind);zp;xframe(1:nwind2)];
    xw = wzp .* xzp;
    Xframe = fft(xw); % FFT
    Xframe=Xframe(select); 
    
    % Calculo del espectrograma en frecuencia midi
    for midi_index=1:muestrasmidi,
        
        kmin=miditobinskmin(midi_index);
        kmax=miditobinskmax(midi_index);
        
        if (kmax-kmin)==0,
            cfreq_amplitudes(midi_index,t) = abs(Xframe(kmin));
        elseif (kmax-kmin)>0,
            cfreq_amplitudes(midi_index,t) = sqrt(sum(abs(Xframe(kmin:kmax)).^2)); % / (kmax-kmin+1));
        end;
    end;
end;

cfreq_amplitudes = cfreq_amplitudes';

% Cierro la puta barra
if barra
    close(h);
end;

return;