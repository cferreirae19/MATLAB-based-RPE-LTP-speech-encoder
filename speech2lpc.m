function [LPC,Ep,RC,LSF,LAR] = speech2lpc(s,p,window,wshift)

% -Funcionalidad: 
%    ·Calcular los coeficientes LPC, RC, LAR y LSF, así como la energía del
%     error de predicción de cada trama.
% -Parámetros de entrada:
%    ·s: Señal de audio
%    ·p: Orden de predicción
%    ·window: Ventana empleada
%    ·wshift: Desplazamiento de la ventana en muestras
% -Parámetros de salida:
%    ·LPC: Matriz con los coeficientes LPC
%    ·Ep: Energía del error de predicción
%    ·RC: Matriz con los coeficientes RC
%    ·LSF: Matriz con los coeficientes LSF
%    ·LAR: Matriz con los coeficientes LAR

for k=1:(length(s)/wshift)  % Iteramos length(s)/wshift veces, es decir, un número de veces igual al número de tramas

    first_sample = (k-1)*wshift + 1;
    last_sample = first_sample + length(window) - 1;

    if (last_sample < length(s))    % Si hay suficientes muestras para crear una trama 

        s_frame = s(first_sample:last_sample);

    else

        s_frame = [s(first_sample:length(s)) zeros(1,(last_sample-length(s)))]; % Si no hay suficientes muestras para crear una trama (por ejemplo, al final) 

    end

    frame = s_frame.*window;    % Enventanamos la trama

    Rcorr = xcorr(frame,frame); % Hallamos todos los coeficientes
    Rcorr = Rcorr(length(frame):length(frame)+p);

    [LPC(k,:),Ep(k),RC(k,:)] = levinson(Rcorr,p);

    LSF(k,:) = poly2lsf(LPC(k,:));

    LAR(k,:) = rc2lar(RC(k,:));

end
