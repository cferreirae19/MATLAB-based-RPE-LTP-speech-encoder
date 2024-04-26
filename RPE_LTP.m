function [xr,eret,ak,bret,Nret] = RPE_LTP(x)

% -Funcionalidad:
%    ·Realizar un codificador para análisis RPE-LTP a partir del trabajo 
%     concatenado de un predictor corto, que eliminará la información del 
%     tracto vocal, y un predictor largo, que eliminará la información 
%     acerca de la periodicidad. En el caso de síntesis RPE-LTP, ambos
%     predictores realizarán los trabajos inversos a sus respectivas etapas
%     en análisis.

num = [1,-0.86];
denom = 1;
s = filter(num,denom,x);  % Filtro de pre-énfasis

z1 = [];    % Condiciones iniciales y finales de filtros
z2 = [];

Lframe = 160;   % Tamaño de trama
wshift = 160;   % Desplazamiento de ventana
window = rectwin(160);  % Ventana
p = 10; % Orden de predicción

Lsubframe = 40; % Tamaño de subtrama

last_subframe = zeros(40,1);  % Última subtrama, para actualizar el vector dp (inicialmente a cero)
last_subframe_r = zeros(40,1);  % Última subtrama recibida, para actualizar el vector drp (inicialmente a cero) 

d = zeros(40,1);    % Señales intermedias del codificador a cero inicialmente
dp = zeros(120,1);
dpN = zeros(40,1);
dpp = zeros(40,1);
e = zeros(40,1);

dr = zeros(40,1);    % Señales intermedias del decodificador a cero inicialmente
drp = zeros(120,1);
drpN = zeros(40,1);
drpp = zeros(40,1);
er = zeros(40,1);

sr_aux = zeros(40,1);
eret=[];
sr = [];
xr = [];

Nret = [];  % Valores que retornaremos de N y b
bret = [];

[ak,~,~,~,~] = speech2lpc(s,p,window,wshift);    % Obtenemos los LPC

for k=1:floor(length(s)/Lsubframe)
    
    first_sample = (k-1)*Lsubframe + 1;

    if(k==floor(length(s)/Lsubframe) && length(s(first_sample:end))<Lsubframe)

        s = [s; zeros(Lsubframe-length(s(first_sample:end)),1)];
    
    end

    % Filtramos cada subtrama con los ak para obtener d
    subframe = s(Lsubframe*(k-1)+1:Lsubframe*k);
    k_4 = ceil(k/4); % Mantenemos los mismos ak (mismo índice) para cada 4 subtramas
    [d,z1] = filter(ak(k_4,:),1,subframe,z1);

    % En cada iteración (para cada subtrama), reorganizamos dp de la siguiente manera:
    %   ·La última subtrama (la más antigua) se perderá (prescindimos de dp(1:40)) 
    %   ·La penúltima subtrama (la "del medio") pasará a ser la más antigua (dp(41:80) pasará  a ser dp(1:40) en nuestro nuevo vector)
    %   ·La antepenúltima subtrama (la más nueva) pasará a ser la "del medio" (dp(81:120) pasará  a ser dp(41:80) en nuestro nuevo vector)
    dp_1_40_new = dp(41:80);
    dp_41_80_new = dp(81:120);

    % Añadimos la última subtrama a dp (pasará a ser dp(81:120) en el nuevo vector dp)
    dp = vertcat(dp_1_40_new,dp_41_80_new,last_subframe);

    % Actualizamos la última subtrama
    last_subframe = d;

    % Hallamos la correlación cruzada entre d y dp, así como su máximo
    R = xcorr(d,dp);
    [Rmax,Rmax_pos] = max(R(40:120));

    % Hallamos el retardo N
    N = Rmax_pos + 39;
    Nret = [Nret,N];    % Almacenamos el valor de N obtenido ahora con los ya existentes

    % Obtenemos la señal dpN (vector de 40 muestras) desplazando dp N muestras
        dpN = dp((121-N):(160-N));


    % Calculamos la ganancia b
    b = Rmax/sum(dpN.^2);
    
    if ~isfinite(b) %ininf(b)||isnan(b)
    b=1;%Eliminamos b si se obtienen valores no coherentes
    end
     
    bret = [bret,b];    % Almacenamos el valor de b obtenido ahora con los ya existentes

    % Amplificamos dpN para obtener dpp
    dpp = b * dpN;

    % Obtenemos e como la resta de d menos dpp
    e = d - dpp;
    
    %Añadimos la subtrama actual al error de predicción
    eret=[eret;e];

    %===DECODIFICACIÓN===%

    % Para decodificación, suponemos er=e
    er = e;

        % En cada iteración (para cada subtrama), reorganizamos drp de la siguiente manera:
    %   ·La última subtrama (la más antigua) se perderá (prescindimos de drp(1:40)) 
    %   ·La penúltima subtrama (la "del medio") pasará a ser la más antigua (drp(41:80) pasará  a ser drp(1:40) en nuestro nuevo vector)
    %   ·La antepenúltima subtrama (la más nueva) pasará a ser la "del medio" (drp(81:120) pasará  a ser drp(41:80) en nuestro nuevo vector)
    drp_1_40_new = drp(41:80);
    drp_41_80_new = drp(81:120);

    % Añadimos la última subtrama a drp (pasará a ser drp(81:120) en el nuevo vector drp)
    drp = vertcat(drp_1_40_new,drp_41_80_new,last_subframe_r);
    

    % Obtenemos la señal drpN (vector de 40 muestras) desplazando drp N muestras
    drpN = drp((121-N):(160-N));

    % Obtenemos drpp amplificando drpN
    drpp = b * drpN;
    % Suma de las señales er y drpp
    
    dr = er + drpp;
    
    % Actualizamos la última subtrama
    last_subframe_r = dr;



    % Filtrado inverso con los ak
    [sr_aux, z2] = filter(1,ak(k_4,:),dr,z2);
    
    fail=d-dr;

    % Añadimos la subtrama actual reconstruida a sr
    sr = [sr;sr_aux];
    

end

% Filtro de de-énfasis
xr = filter(denom,num,sr);

end