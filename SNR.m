
function S = SNR(x,xq)

% -Funcionalidad:
%    ·Calcular la relación señal a ruido (SNR) de una señal dada
% -Parámetros de entrada:
%    ·x: Señal original
%    ·xq: Señal cuantificada
% -Parámetros de salida:
%    ·S: Relación señal a ruido

S = 10 * log10(sum(x.^2)/sum((x-xq).^2));

end

