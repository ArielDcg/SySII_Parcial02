%% TEST ONE FREQ - PRUEBA DE UNA FRECUENCIA INDIVIDUAL
% Script para probar el ecualizador con una sola frecuencia a la vez
% Muestra las respuestas en magnitud y fase de todos los filtros
% y los espectros de entrada y salida de cada filtro

clear all;
close all;
clc;

%% Parámetros
fs = 64000; % Hz

%% Seleccionar frecuencia a probar
fprintf('Seleccione la frecuencia a probar:\n');
fprintf('  1. Graves       (102.4 Hz)\n');
fprintf('  2. Medios bajos (512 Hz)\n');
fprintf('  3. Medios altos (2560 Hz)\n');
fprintf('  4. Agudos       (12800 Hz)\n\n');

opcion = input('Opción [1-4]: ');

switch opcion
    case 1
        fc_test = 102.4;
        nombre_banda = 'Graves';
    case 2
        fc_test = 512;
        nombre_banda = 'Medios Bajos';
    case 3
        fc_test = 2560;
        nombre_banda = 'Medios Altos';
    case 4
        fc_test = 12800;
        nombre_banda = 'Agudos';
    otherwise
        error('Opción inválida');
end

fprintf('\nProbando banda: %s (fc = %.1f Hz)\n\n', nombre_banda, fc_test);

%% Diseñar los 4 filtros
fprintf('Diseñando los 4 filtros del ecualizador...\n');

[b1, a1, fc1, tipo1, N1, H_freq1] = filtro_graves();
[b2, a2, fc2, tipo2, N2, H_freq2] = filtro_medios_bajos();
[b3, a3, fc3, tipo3, N3, H_freq3] = filtro_medios_altos();
[b4, a4, fc4, tipo4, N4, H_freq4] = filtro_agudos();

filtros_b = {b1, b2, b3, b4};
filtros_a = {a1, a2, a3, a4};
filtros_fc = [fc1, fc2, fc3, fc4];
filtros_tipo = {tipo1, tipo2, tipo3, tipo4};
filtros_N = [N1, N2, N3, N4];
H_freqs = {H_freq1, H_freq2, H_freq3, H_freq4};

fprintf('Filtros diseñados exitosamente.\n\n');

%% Generar señal de prueba (tono puro en la frecuencia seleccionada)
duracion = 1; % segundo
t = 0:1/fs:duracion-1/fs;
audio_test = sin(2*pi*fc_test*t);

% Normalizar a rango de 16 bits
audio_test = audio_test / max(abs(audio_test));
audio_test = audio_test * (2^15 - 1);

fprintf('Señal de prueba: tono puro a %.1f Hz\n', fc_test);
fprintf('Duración: %.1f s, Muestras: %d\n\n', duracion, length(audio_test));

%% Calcular espectro de entrada
N_fft = 2^nextpow2(length(audio_test));
X_in = fft(audio_test, N_fft);
f_fft = (0:N_fft-1) * fs / N_fft;
idx_pos = 1:floor(N_fft/2)+1;
f_espectro = f_fft(idx_pos);
mag_in_dB = 20*log10(abs(X_in(idx_pos)) + eps);

%% GRÁFICA 1: Respuestas en magnitud y fase de todos los filtros
fprintf('Generando gráfica 1: Respuestas en magnitud y fase...\n');

figure('Position', [50, 50, 1400, 800]);

% Magnitud
subplot(2,1,1);
colores = ['b', 'r', 'g', 'm'];
hold on;
for i = 1:4
    semilogx(H_freqs{i}.f, H_freqs{i}.mag_dB, colores(i), 'LineWidth', 2, ...
             'DisplayName', sprintf('Banda %d: %s (%.1f Hz)', i, filtros_tipo{i}, filtros_fc(i)));
end
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Respuestas en Magnitud de Todos los Filtros');
legend('Location', 'southwest');
xlim([10 fs/2]);
ylim([-80 5]);
yline(-3, '--k', 'LineWidth', 1);
yline(-15, '--k', 'LineWidth', 0.5);

% Marcar frecuencia de prueba
xline(fc_test, '--r', sprintf('%.1f Hz', fc_test), 'LineWidth', 1.5);

% Fase
subplot(2,1,2);
hold on;
for i = 1:4
    semilogx(H_freqs{i}.f, H_freqs{i}.fase_deg, colores(i), 'LineWidth', 2, ...
             'DisplayName', sprintf('Banda %d', i));
end
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
title('Respuestas en Fase de Todos los Filtros');
legend('Location', 'best');
xlim([10 fs/2]);

% Marcar frecuencia de prueba
xline(fc_test, '--r', sprintf('%.1f Hz', fc_test), 'LineWidth', 1.5);

sgtitle(sprintf('Respuestas de Frecuencia - Banda de Prueba: %s', nombre_banda), ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, sprintf('respuestas_filtros_%s.png', strrep(nombre_banda, ' ', '_')));
fprintf('Guardado: respuestas_filtros_%s.png\n\n', strrep(nombre_banda, ' ', '_'));

%% GRÁFICA 2: Espectros de entrada y salida de cada filtro
fprintf('Generando gráfica 2: Espectros de entrada y salida...\n');

figure('Position', [100, 100, 1400, 900]);

for i = 1:4
    subplot(2, 2, i);

    % Espectro de entrada (mismo para todos)
    semilogx(f_espectro, mag_in_dB, 'b', 'LineWidth', 1.5, 'DisplayName', 'Entrada');
    hold on;

    % Espectro de salida (por ahora igual a entrada, pendiente procesamiento)
    % TODO: Calcular salida real cuando se implemente procesamiento de audio
    mag_out_dB = mag_in_dB; % Temporal
    semilogx(f_espectro, mag_out_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Salida');

    % Superponer respuesta del filtro
    semilogx(H_freqs{i}.f, H_freqs{i}.mag_dB, colores(i), 'LineWidth', 1, ...
             'DisplayName', 'H(jω)', 'LineStyle', ':');

    grid on;
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Banda %d: %s (fc=%.1f Hz)', i, filtros_tipo{i}, filtros_fc(i)));
    legend('Location', 'best');
    xlim([10 fs/2]);
    ylim([-100 100]);

    % Marcar frecuencia de prueba
    xline(fc_test, '--k', 'LineWidth', 1);
end

sgtitle(sprintf('Espectros de Entrada y Salida - Señal de Prueba: %.1f Hz', fc_test), ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, sprintf('espectros_filtros_%s.png', strrep(nombre_banda, ' ', '_')));
fprintf('Guardado: espectros_filtros_%s.png\n\n', strrep(nombre_banda, ' ', '_'));

%% Resumen
fprintf('=========================================================================\n');
fprintf('                         RESUMEN DE LA PRUEBA\n');
fprintf('=========================================================================\n\n');

fprintf('Banda probada: %s (fc = %.1f Hz)\n', nombre_banda, fc_test);
fprintf('Señal de prueba: Tono puro a %.1f Hz\n\n', fc_test);

fprintf('Filtros del ecualizador:\n');
for i = 1:4
    fprintf('  Banda %d: %s (fc=%.1f Hz, Orden=%d)\n', ...
            i, filtros_tipo{i}, filtros_fc(i), filtros_N(i));
end
fprintf('\n');

fprintf('Gráficas generadas:\n');
fprintf('  1. Respuestas en magnitud y fase de todos los filtros\n');
fprintf('  2. Espectros de entrada y salida de cada filtro\n\n');

fprintf('NOTA: Procesamiento de audio pendiente de implementación.\n');
fprintf('      Las salidas mostradas son temporales.\n\n');

fprintf('=========================================================================\n');
fprintf('                      PRUEBA COMPLETADA\n');
fprintf('=========================================================================\n\n');
