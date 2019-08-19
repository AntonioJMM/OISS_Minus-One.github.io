function [yo,fo,to] = sg_slow(x,nfft,Fs,win,noverlap,barra)
% S = sg_slow(B,NFFT,Fs,win,NOVERLAP)
% win must be created like win = hanning(Nsamples)
% All parameters like SPECGRAM

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

y=zeros(length(select),nframes);

% X = zeros(length(select),nframes);
for t=1:nframes,
    % Actualizo la barra
    if barra && mod(t,10)==0,
        waitbar(t/nframes,h,['Trama ' num2str(t) ' de ' num2str(nframes)]);
    end;
    
    xframe = x(frame_index(t):frame_index(t)+nwind-1); % extract frame of input data
    xzp = [xframe(nwind2+1:nwind);zp;xframe(1:nwind2)];
    xw = wzp .* xzp;
    Xframe = fft(xw); % FFT
    Xframe=Xframe(select); 
    
    y(:,t)=Xframe;
end;

% Cierro la puta barra
if barra
    close(h);
end;

% f = (select - 1)'*Fs/nfft;
% colindex = 1 + (0:(nframes-1))*(nwind-noverlap);
% t = (colindex-1)'/Fs;
% if nargout == 1,
%     yo = y;
% elseif nargout == 2,
%     yo = y;
%     fo = f;
% elseif nargout == 3,
%     yo = y;
%     fo = f;
%     to = t;
% end
yo = y;

return;
