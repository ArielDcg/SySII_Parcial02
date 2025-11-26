function [b, a, fc, tipo_filtro, N, H_freq] = filtro_graves()
%FILTRO_GRAVES Diseña filtro pasabajos Butterworth para banda de graves (102.4 Hz)
%   Salidas:
%       b, a: Coeficientes del filtro H(s)
%       fc: Frecuencia central (102.4 Hz)
%       tipo_filtro: 'PASABAJOS'
%       N: Orden del filtro
%       H_freq: Estructura con respuesta en frecuencia

    %% Especificaciones
    fc = 102.4; % Hz
    f_central = [102.4, 512, 2560, 12800];
    fs = 64000; % Hz
    Rp = 3;   % Atenuación en banda de paso (dB)
    As = 15;  % Atenuación en banda de rechazo (dB)
    tipo_filtro = 'PASABAJOS';

    %% Frecuencias límite
    fp = sqrt(f_central(1) * f_central(2));  % -3dB
    fs_rechazo = f_central(2);               % -15dB
    Wp = 2*pi*fp;
    Ws = 2*pi*fs_rechazo;

    %% Calcular orden del filtro
    [N, Wc] = buttord(Wp, Ws, Rp, As, 's');

    %% Diseñar filtro Butterworth
    [b, a] = butter(N, Wc, 's');

    %% Calcular respuesta en frecuencia
    w_plot = logspace(log10(2*pi*10), log10(2*pi*fs/2), 2000);
    [H, w_resp] = freqs(b, a, w_plot);
    f_plot = w_resp/(2*pi);
    mag_dB = 20*log10(abs(H));
    fase_deg = angle(H)*180/pi;

    % Guardar en estructura
    H_freq.w = w_plot;
    H_freq.f = f_plot;
    H_freq.H = H;
    H_freq.mag_dB = mag_dB;
    H_freq.fase_deg = fase_deg;

end
