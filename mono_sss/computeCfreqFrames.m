function [X_ft]=computeCfreqFrames(x,fft_params,active_frames,draw)
% [X_ft,miditobins]=computeCfreqSeq(x,fft_params,draw)
% Define the default fft_paramseters from the input signal
%
% Inputs:
% x  - signal
% fft_params - input parameters 
% draw - show process bar (default 0) 
%
% Outputs:
% X_ft - amplitude "MIDI" spectrogram (mnotesXframes)
% miditobins - Relation MIDI note & MIDI bins (2Xmnotes)
%
% Julio Carabias y Francisco Rodriguez. Fall 2012

if nargin<3,
    draw = 0;
end

% Inicializaciones
fftsize = fft_params.fftsize;
nframes = numel(active_frames);

% Genero escala MIDI
if isfield(fft_params,'miditobins') && isfield(fft_params,'muestrasmidi')
    miditobins = fft_params.miditobins;
	muestrasmidi = fft_params.muestrasmidi;
else
    [miditobins,muestrasmidi]=computeCfreqInit(fft_params);
end

% Inicializacion de variables de salida
X_ft=zeros(numel(active_frames),muestrasmidi);
X_ft = X_ft';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculo el spectrograma secuencialmente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowsize = fft_params.windowsize;
hopsize = fft_params.hopsize;
nfft=2*fftsize;
win = hanning(windowsize,'periodic');  % Quitado SQRT a la ventana. FCO. 24/11/11
noverlap = windowsize-hopsize;

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

%nframes = nframes;
nframes = fix((nx-noverlap)/(nwind-noverlap));
frame_index = 1 + (0:(nframes-1))*(nwind-noverlap);
if length(x)<(nwind+frame_index(nframes)-1)
    x(end+1:nwind+frame_index(nframes)-1) = 0;   % zero-pad x
end

% La puta draw
if draw,
    h = waitbar(0,'Please wait...','Name','Computing Spectrogram ...');
end;

% X = zeros(length(select),nframes);
for t_iter= 1:numel(active_frames),
    t = active_frames(t_iter);
    
    % Actualizo la draw
    if draw && mod(t,10)==0,
        waitbar(t_iter/nframes,h,['Frame ' num2str(t_iter) ' of ' num2str(nframes)]);
    end;
    
    xframe = x(frame_index(t):frame_index(t)+nwind-1); % extract frame of input data
    xzp = [xframe(nwind2+1:nwind);zp;xframe(1:nwind2)];
    xw = wzp .* xzp;
    Xframe = fft(xw); % FFT
    Xframe=Xframe(select); 
    
    % Calculo del espectrograma en frecuencia midi
    for midi_index=1:muestrasmidi,
        
        kmin=miditobins(1,midi_index);
        kmax=miditobins(2,midi_index);
        
        if (kmax-kmin)==0,
            X_ft(midi_index,t) = abs(Xframe(kmin));
        elseif (kmax-kmin)>0,
            X_ft(midi_index,t) = sqrt(sum(abs(Xframe(kmin:kmax)).^2)); % / (kmax-kmin+1));
        end;
    end;
end;

% Cierro la puta draw
if draw
    close(h);
end;

return;