function [b, a, fc, tipo_filtro, N, H_freq] = filtro_agudos()
%FILTRO_AGUDOS Diseña filtro pasaaltas Butterworth para banda de agudos (12800 Hz)
%   Salidas:
%       b, a: Coeficientes del filtro H(s)
%       fc: Frecuencia central (12800 Hz)
%       tipo_filtro: 'PASAALTAS'
%       N: Orden del filtro
%       H_freq: Estructura con respuesta en frecuencia

    %% Especificaciones
    fc = 12800; % Hz
    f_central = [102.4, 512, 2560, 12800];
    fs = 64000; % Hz
    Rp = 3;   % Atenuación en banda de paso (dB)
    As = 15;  % Atenuación en banda de rechazo (dB)
    tipo_filtro = 'PASAALTAS';

    %% Frecuencias límite
    fp = sqrt(f_central(3) * f_central(4));  % -3dB
    fs_rechazo = f_central(3);               % -15dB
    Wp = 2*pi*fp;
    Ws = 2*pi*fs_rechazo;

    %% Calcular orden del filtro
    [N, Wc] = buttord(Wp, Ws, Rp, As, 's');

    %% Diseñar filtro Butterworth pasaaltas
    [b_lp, a_lp] = butter(N, 1, 's');
    [b, a] = lp2hp(b_lp, a_lp, Wc);

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
