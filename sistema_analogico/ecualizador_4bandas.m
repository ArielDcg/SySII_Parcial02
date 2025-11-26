function [audio_out, filtros, espectros] = ecualizador_4bandas(audio_in, ganancias, fs)
%ECUALIZADOR_4BANDAS Ecualiza audio usando 4 filtros analógicos Butterworth
%   Entradas:
%       audio_in: Señal de audio de entrada (16 bits, mono)
%       ganancias: Vector de 4 ganancias [G1, G2, G3, G4] para cada banda
%       fs: Frecuencia de muestreo (default: 64000 Hz)
%   Salidas:
%       audio_out: Señal de audio ecualizada
%       filtros: Estructura con información de los 4 filtros
%       espectros: Estructura con espectros de entrada y salida

    %% Validar entradas
    if nargin < 3
        fs = 64000; % Hz por defecto
    end

    if nargin < 2
        error('Se requieren al menos 2 argumentos: audio_in y ganancias');
    end

    if length(ganancias) ~= 4
        error('El vector de ganancias debe tener 4 elementos [G1, G2, G3, G4]');
    end

    % Asegurar que audio_in es monofonico (columna)
    if size(audio_in, 2) > 1
        warning('Audio de entrada tiene múltiples canales. Convirtiendo a mono...');
        audio_in = mean(audio_in, 2);
    end

    % Convertir a columna si es necesario
    if size(audio_in, 1) == 1
        audio_in = audio_in(:);
    end

    %% Diseñar los 4 filtros
    fprintf('Diseñando filtros del ecualizador...\n');

    % Graves (102.4 Hz)
    [b1, a1, fc1, tipo1, N1, H_freq1] = filtro_graves();

    % Medios bajos (512 Hz)
    [b2, a2, fc2, tipo2, N2, H_freq2] = filtro_medios_bajos();

    % Medios altos (2560 Hz)
    [b3, a3, fc3, tipo3, N3, H_freq3] = filtro_medios_altos();

    % Agudos (12800 Hz)
    [b4, a4, fc4, tipo4, N4, H_freq4] = filtro_agudos();

    % Guardar información de filtros
    filtros.b = {b1, b2, b3, b4};
    filtros.a = {a1, a2, a3, a4};
    filtros.fc = [fc1, fc2, fc3, fc4];
    filtros.tipo = {tipo1, tipo2, tipo3, tipo4};
    filtros.N = [N1, N2, N3, N4];
    filtros.H_freq = {H_freq1, H_freq2, H_freq3, H_freq4};

    fprintf('Filtros diseñados:\n');
    fprintf('  Banda 1: %s (fc=%.1f Hz, N=%d)\n', tipo1, fc1, N1);
    fprintf('  Banda 2: %s (fc=%.0f Hz, N=%d)\n', tipo2, fc2, N2);
    fprintf('  Banda 3: %s (fc=%.0f Hz, N=%d)\n', tipo3, fc3, N3);
    fprintf('  Banda 4: %s (fc=%.0f Hz, N=%d)\n', tipo4, fc4, N4);

    %% Calcular espectro de entrada
    N_fft = 2^nextpow2(length(audio_in));
    X_in = fft(audio_in, N_fft);
    f_fft = (0:N_fft-1) * fs / N_fft;

    % Solo frecuencias positivas
    idx_pos = 1:floor(N_fft/2)+1;
    espectros.f = f_fft(idx_pos);
    espectros.X_in = X_in(idx_pos);
    espectros.mag_in_dB = 20*log10(abs(X_in(idx_pos)) + eps);

    %% Procesar audio (preparado para 16 bits mono)
    % Por ahora, inicializar salida (procesamiento de audio pendiente)
    audio_out = zeros(size(audio_in));

    % TODO: Implementar procesamiento de audio con filtros analógicos
    % Esto requerirá conversión de filtros analógicos a digitales
    % usando bilinear, impulse invariance, u otro método

    fprintf('\nNOTA: Procesamiento de audio aún no implementado.\n');
    fprintf('      Los filtros H(s) están diseñados y listos.\n');
    fprintf('      Se requiere conversión a filtros digitales H(z).\n\n');

    %% Calcular espectro de salida (temporal - será reemplazado)
    X_out = fft(audio_out, N_fft);
    espectros.X_out = X_out(idx_pos);
    espectros.mag_out_dB = 20*log10(abs(X_out(idx_pos)) + eps);

    % Guardar ganancias usadas
    espectros.ganancias = ganancias;

end
