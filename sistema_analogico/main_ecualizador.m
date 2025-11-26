%% MAIN ECUALIZADOR - SISTEMA ANALÓGICO 4 BANDAS
% Script principal que prueba el ecualizador con diferentes ganancias
% y genera espectros de magnitud de entrada y salida

clear all;
close all;
clc;

%% Parámetros
fs = 64000; % Hz
duracion = 1; % segundo
t = 0:1/fs:duracion-1/fs;

%% Generar señal de prueba
% Señal de prueba: suma de tonos en diferentes bandas
f_test = [100, 500, 2500, 12000]; % Hz (una frecuencia por banda)
audio_test = zeros(size(t));
for i = 1:length(f_test)
    audio_test = audio_test + sin(2*pi*f_test(i)*t);
end

% Normalizar a rango de 16 bits
audio_test = audio_test / max(abs(audio_test));
audio_test = audio_test * (2^15 - 1); % Escalar a 16 bits

fprintf('Señal de prueba generada:\n');
fprintf('  Duración: %.1f s\n', duracion);
fprintf('  Frecuencias: [%d, %d, %d, %d] Hz\n', f_test(1), f_test(2), f_test(3), f_test(4));
fprintf('  Muestras: %d\n\n', length(audio_test));

%% Vector de ganancias a probar
ganancias_vector = [0.25, 0.5, 1, 2];
n_ganancias = length(ganancias_vector);

%% Probar ecualizador con diferentes ganancias (todas las bandas con la misma ganancia)
fprintf('=========================================================================\n');
fprintf('          PRUEBA DE ECUALIZADOR CON DIFERENTES GANANCIAS\n');
fprintf('=========================================================================\n\n');

resultados = cell(n_ganancias, 1);

for i = 1:n_ganancias
    G = ganancias_vector(i);
    fprintf('Prueba %d: Ganancia = %.2f (todas las bandas)\n', i, G);
    fprintf('-------------------------------------------------------------------------\n');

    % Todas las bandas con la misma ganancia
    ganancias = [G, G, G, G];

    % Ejecutar ecualizador
    [audio_out, filtros, espectros] = ecualizador_4bandas(audio_test, ganancias, fs);

    % Guardar resultados
    resultados{i}.ganancias = ganancias;
    resultados{i}.audio_out = audio_out;
    resultados{i}.filtros = filtros;
    resultados{i}.espectros = espectros;

    fprintf('\n');
end

%% Generar gráficas de espectros de magnitud
fprintf('Generando gráficas de espectros de magnitud...\n\n');

figure('Position', [50, 50, 1400, 900]);

for i = 1:n_ganancias
    G = ganancias_vector(i);

    % Subplot para cada ganancia
    subplot(2, 2, i);

    % Espectro de entrada
    semilogx(resultados{i}.espectros.f, resultados{i}.espectros.mag_in_dB, ...
             'b', 'LineWidth', 1.5, 'DisplayName', 'Entrada');
    hold on;

    % Espectro de salida
    semilogx(resultados{i}.espectros.f, resultados{i}.espectros.mag_out_dB, ...
             'r', 'LineWidth', 1.5, 'DisplayName', 'Salida');

    % Marcar frecuencias de prueba
    for j = 1:length(f_test)
        xline(f_test(j), '--k', 'LineWidth', 0.5);
    end

    grid on;
    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Ganancia = %.2f', G));
    legend('Location', 'best');
    xlim([10 fs/2]);
    ylim([-100 100]);
end

sgtitle('Espectros de Magnitud - Entrada y Salida del Sistema', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Guardar figura
saveas(gcf, 'espectros_magnitud_ganancias.png');
fprintf('Figura guardada: espectros_magnitud_ganancias.png\n\n');

%% Resumen
fprintf('=========================================================================\n');
fprintf('                           RESUMEN DE PRUEBAS\n');
fprintf('=========================================================================\n\n');

fprintf('Ganancias probadas: [%.2f, %.2f, %.2f, %.2f]\n', ganancias_vector);
fprintf('Filtros diseñados:\n');
fprintf('  Banda 1 - Graves: fc=%.1f Hz (%s, N=%d)\n', ...
        resultados{1}.filtros.fc(1), resultados{1}.filtros.tipo{1}, resultados{1}.filtros.N(1));
fprintf('  Banda 2 - Medios bajos: fc=%.0f Hz (%s, N=%d)\n', ...
        resultados{1}.filtros.fc(2), resultados{1}.filtros.tipo{2}, resultados{1}.filtros.N(2));
fprintf('  Banda 3 - Medios altos: fc=%.0f Hz (%s, N=%d)\n', ...
        resultados{1}.filtros.fc(3), resultados{1}.filtros.tipo{3}, resultados{1}.filtros.N(3));
fprintf('  Banda 4 - Agudos: fc=%.0f Hz (%s, N=%d)\n', ...
        resultados{1}.filtros.fc(4), resultados{1}.filtros.tipo{4}, resultados{1}.filtros.N(4));
fprintf('\n');

fprintf('=========================================================================\n');
fprintf('                      PROCESAMIENTO COMPLETADO\n');
fprintf('=========================================================================\n\n');
