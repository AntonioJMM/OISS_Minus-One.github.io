function ventanas=nmf_generateWindExcitINC(param)

ventanas=zeros((param.midi_max-param.midi_min+1)*param.midi_inc ,param.muestrasflog);

pitches=param.midi_min:1/param.midi_inc:param.midi_max+(param.midi_inc-1)*(1/param.midi_inc);
frecuencias=((2.^((pitches-69)/12))*440);
 
%ventana = hanning(param.windowsize,'periodic');
ventana = sqrt(hanning(param.windowsize,'periodic'));  % SQRT quitado a la ventana. FCO. 24/11/11

t=0:1/param.fs:param.t_trama-(1/param.fs);
%x = sin(2*pi*frecuencias'*t).*repmat(ventana,size(frecuencias),1)';
x = bsxfun(@times,sin(2*pi*frecuencias'*t),ventana');

n_tramas=fix((size(x,2)-param.windowsize)/param.salto)+1;
param.n_tramas= n_tramas;

for k=1:size(ventanas,1),
        %X = sg(x(k,:),2*param.longitud_frec,param.fs,ventana,param.windowsize-param.salto);
        %[cfreq_amplitudes,miditobins] = modelo_rapido_suma(X,param.fs,param.longitud_frec,n_tramas,param.midi_inc); %Modelo rapido suma de intervalos
        cfreq_amplitudes = computeCfreqFast(x(k,:),param);
        cfreq_amplitudes = cfreq_amplitudes';
        ventanas(k,:) = cfreq_amplitudes(:,round(size(cfreq_amplitudes,2)/2))';
        ventanas(k,:) = ventanas(k,:)./max(ventanas(k,:));
end;

return;
