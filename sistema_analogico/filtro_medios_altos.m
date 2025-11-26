function [b, a, fc, tipo_filtro, N, H_freq] = filtro_medios_altos()
%FILTRO_MEDIOS_ALTOS Diseña filtro pasabanda Butterworth para medios altos (2560 Hz)
%   Salidas:
%       b, a: Coeficientes del filtro H(s)
%       fc: Frecuencia central (2560 Hz)
%       tipo_filtro: 'PASABANDA'
%       N: Orden del filtro
%       H_freq: Estructura con respuesta en frecuencia

    %% Especificaciones
    fc = 2560; % Hz
    f_central = [102.4, 512, 2560, 12800];
    idx_banda = 3;
    fs = 64000; % Hz
    Rp = 3;   % Atenuación en banda de paso (dB)
    As = 15;  % Atenuación en banda de rechazo (dB)
    tipo_filtro = 'PASABANDA';

    %% Frecuencias límite
    f1_3dB = sqrt(f_central(idx_banda-1) * f_central(idx_banda));
    f2_3dB = sqrt(f_central(idx_banda) * f_central(idx_banda+1));

    w1 = 2*pi*f1_3dB;
    w2 = 2*pi*f2_3dB;
    ws1 = 2*pi*f_central(idx_banda-1);
    ws2 = 2*pi*f_central(idx_banda+1);

    Wo = sqrt(w1*w2);
    BW = w2 - w1;

    %% Calcular orden del filtro
    [N_inferior, ~] = buttord(w1, ws1, Rp, As, 's');
    [N_superior, ~] = buttord(w2, ws2, Rp, As, 's');
    N = max(N_inferior, N_superior);

    %% Diseñar filtro Butterworth pasabanda
    [b_lp, a_lp] = butter(N, 1, 's');
    [b, a] = lp2bp(b_lp, a_lp, Wo, BW);

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
