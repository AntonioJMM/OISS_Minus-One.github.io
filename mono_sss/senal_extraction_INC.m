function [senales]=senal_extraction_INC(x,param,A_ptj,S_fpj)
% [senales]=senal_extraction_INC(x,param,channel,A_ptj,S_fpj,M_js)
% Perform the signal source separation using Wiener Filtering.
%
% Input arguments:
%	x --> input signal
%   params --> FFT parameters
%   channel --> Selected channel
%   A_ptj --> NMF Gains
%   S_fpj --> NMF Basis Functions
%
% Output:
%   senales --> separated signals

% Inicializaciones
j_max = size(A_ptj,3);

% Definicion de la salida
senales = zeros(length(x),j_max);
Xe_ftj = tprod(S_fpj,[1 -1 3],A_ptj,[-1 2 3]);

%%% La separacion por casilla trama-midi
denominador = (sum(Xe_ftj.^2,3));
Xe_ftj_relativo = (Xe_ftj.^2) ./ repmat(denominador+realmin,[1,1,j_max]);
Xe_ftj_relativo(repmat(denominador==0,1,j_max))=0;

% Valores de configuracion
ventana = sqrt(hanning(param.windowsize,'periodic'));

% Espectrogramas
X=sg(x(:,1),2*param.longitud_frec,param.fs,ventana,param.windowsize-param.salto);
%X=sg_slow(x(:,1),2*param.longitud_frec,param.fs,ventana,param.windowsize-param.salto,1);
X_r = zeros(size(X,1),size(X,2),j_max);

% Calculo del espectrograma en frecuencia midi INCREMENTADA a 1/4 midi
for inc_midi=param.flog_first:param.muestrasflog
    kmin=param.flog2bins(1,inc_midi);
    kmax=param.flog2bins(2,inc_midi);
    
    if (kmax-kmin)>=0
        for jj=1:j_max
            X_r(kmin:kmax,1:size(Xe_ftj,2),jj) = X(kmin:kmax,1:size(Xe_ftj,2)) .* repmat(Xe_ftj_relativo(inc_midi,:,jj),kmax-kmin+1,1);
        end
    end
end

% Para cada instrumento
for jj=1:j_max
    % Espectrograma inverso
    senal_act = invspecgram(X_r(:,:,jj),2*param.longitud_frec,param.fs,ventana(:),param.windowsize-param.salto); % Ventana vector columna
    senales(1:length(senal_act),jj) = senal_act;
end

return