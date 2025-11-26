%% ECUALIZADOR DE AUDIO ANALÓGICO - 4 BANDAS
% Diseña los 4 filtros analógicos Butterworth para el ecualizador

clear all;
close all;
clc;

%% Parámetros del sistema
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('         DISEÑO COMPLETO DEL ECUALIZADOR ANALÓGICO - 4 BANDAS\n');
fprintf('=========================================================================\n\n');

% Frecuencias centrales de las 4 bandas
f_central = [102.4, 512, 2560, 12800]; % Hz
fs = 64000; % Hz
n_bandas = length(f_central);

fprintf('Frecuencias centrales: [%.1f, %.0f, %.0f, %.0f] Hz\n', ...
        f_central(1), f_central(2), f_central(3), f_central(4));
fprintf('Frecuencia de muestreo: %d Hz\n\n', fs);

%% Diseño de los 4 filtros
fprintf('-------------------------------------------------------------------------\n');
fprintf('                  DISEÑO DE LOS 4 FILTROS DEL ECUALIZADOR\n');
fprintf('-------------------------------------------------------------------------\n\n');

% Arrays para guardar los coeficientes de cada filtro
coeficientes_b = cell(n_bandas, 1);
coeficientes_a = cell(n_bandas, 1);
ordenes = zeros(n_bandas, 1);

% Diseñar cada filtro
for i = 1:n_bandas
    fprintf('*** BANDA %d: fc = %.1f Hz ***\n\n', i, f_central(i));

    % Llamar a la función de diseño (sin gráficas para no saturar)
    disenar_filtro_analogico(f_central(i), false);

    % Cargar los coeficientes guardados
    if f_central(i) == 102.4
        nombre_mat = 'filtro_butterworth_102p4Hz.mat';
    else
        nombre_mat = sprintf('filtro_butterworth_%.0fHz.mat', f_central(i));
    end

    datos = load(nombre_mat);
    coeficientes_b{i} = datos.b;
    coeficientes_a{i} = datos.a;
    ordenes(i) = length(datos.a) - 1;

    fprintf('-------------------------------------------------------------------------\n\n');
end

%% Resumen del diseño
fprintf('=========================================================================\n');
fprintf('                      RESUMEN DEL ECUALIZADOR\n');
fprintf('=========================================================================\n\n');

fprintf('FILTROS DISEÑADOS:\n');
fprintf('  Banda 1: fc = %.1f Hz  -> PASABAJOS  (Orden: %d)\n', f_central(1), ordenes(1));
fprintf('  Banda 2: fc = %.0f Hz   -> PASABANDA  (Orden: %d)\n', f_central(2), ordenes(2));
fprintf('  Banda 3: fc = %.0f Hz  -> PASABANDA  (Orden: %d)\n', f_central(3), ordenes(3));
fprintf('  Banda 4: fc = %.0f Hz -> PASAALTAS  (Orden: %d)\n\n', f_central(4), ordenes(4));

%% Gráfica comparativa de todos los filtros
fprintf('Generando gráfica comparativa de respuestas en frecuencia...\n\n');

% Vector de frecuencias
w_plot = logspace(log10(2*pi*10), log10(2*pi*fs/2), 2000);
f_plot = w_plot/(2*pi);

% Calcular respuestas en frecuencia
H_filtros = zeros(n_bandas, length(w_plot));
for i = 1:n_bandas
    [H, ~] = freqs(coeficientes_b{i}, coeficientes_a{i}, w_plot);
    H_filtros(i, :) = 20*log10(abs(H));
end

% Crear figura comparativa
figure('Position', [50, 50, 1400, 900]);

% Subplot 1: Todas las respuestas individuales
subplot(2,2,1);
colores = ['b', 'r', 'g', 'm'];
hold on;
for i = 1:n_bandas
    semilogx(f_plot, H_filtros(i, :), colores(i), 'LineWidth', 2, ...
             'DisplayName', sprintf('Banda %d (%.1f Hz)', i, f_central(i)));
end
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Respuestas Individuales de los 4 Filtros');
legend('Location', 'southwest');
xlim([10 fs/2]);
ylim([-80 5]);
yline(-3, '--k', 'LineWidth', 1);
yline(-15, '--k', 'LineWidth', 0.5);

% Subplot 2: Suma de todas las respuestas (ecualización plana G=1)
subplot(2,2,2);
H_suma = sum(10.^(H_filtros/20), 1); % Suma en lineal
H_suma_dB = 20*log10(H_suma);
semilogx(f_plot, H_suma_dB, 'k', 'LineWidth', 2);
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Respuesta Total del Ecualizador (G1=G2=G3=G4=1)');
xlim([10 fs/2]);
ylim([-10 10]);
yline(0, '--r', 'LineWidth', 1.5);

% Subplot 3: Vista detallada de cada filtro alrededor de su fc
subplot(2,2,3);
hold on;
for i = 1:n_bandas
    semilogx(f_plot, H_filtros(i, :), colores(i), 'LineWidth', 2);
    % Marcar frecuencia central
    plot(f_central(i), 0, 'o', 'Color', colores(i), 'MarkerSize', 8, 'LineWidth', 2);
end
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Detalle: Frecuencias Centrales');
xlim([50 20000]);
ylim([-40 5]);
yline(-3, '--k', 'LineWidth', 1);

% Subplot 4: Cobertura de frecuencias (vista lineal de magnitud)
subplot(2,2,4);
hold on;
for i = 1:n_bandas
    H_lin = 10.^(H_filtros(i, :)/20);
    semilogx(f_plot, H_lin, colores(i), 'LineWidth', 2, ...
             'DisplayName', sprintf('Banda %d', i));
end
grid on;
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (lineal)');
title('Respuestas en Escala Lineal');
legend('Location', 'best');
xlim([10 fs/2]);
ylim([0 1.2]);
yline(1, '--k', 'LineWidth', 1);
yline(0.707, '--r', 'LineWidth', 1); % -3dB en lineal

sgtitle('Ecualizador de Audio Analógico - 4 Bandas', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Guardar resultados
saveas(gcf, 'ecualizador_analogico_completo.png');
fprintf('Gráfica guardada: ecualizador_analogico_completo.png\n');

% Guardar todos los coeficientes en un solo archivo
save('ecualizador_analogico.mat', 'coeficientes_b', 'coeficientes_a', ...
     'ordenes', 'f_central', 'fs');
fprintf('Coeficientes guardados: ecualizador_analogico.mat\n\n');

%% Análisis de cobertura de frecuencias
fprintf('=========================================================================\n');
fprintf('                  ANÁLISIS DE COBERTURA DE FRECUENCIAS\n');
fprintf('=========================================================================\n\n');

% Puntos de cruce entre filtros adyacentes (-3dB)
fprintf('PUNTOS DE CRUCE ENTRE FILTROS (a -3dB):\n');
for i = 1:n_bandas-1
    fc_cruce = sqrt(f_central(i) * f_central(i+1));
    fprintf('  Entre Banda %d y %d: %.2f Hz\n', i, i+1, fc_cruce);
end
fprintf('\n');

% Verificar planicidad de la respuesta total
fprintf('PLANICIDAD DE LA RESPUESTA TOTAL (G=1 para todas las bandas):\n');
banda_audio = (f_plot >= 20) & (f_plot <= 20000);
H_suma_banda_audio = H_suma_dB(banda_audio);
max_ripple = max(H_suma_banda_audio) - min(H_suma_banda_audio);
fprintf('  Rango de frecuencias de audio: 20 Hz - 20 kHz\n');
fprintf('  Magnitud máxima: %.2f dB\n', max(H_suma_banda_audio));
fprintf('  Magnitud mínima: %.2f dB\n', min(H_suma_banda_audio));
fprintf('  Ripple (variación): %.2f dB\n\n', max_ripple);

if max_ripple <= 3
    fprintf('  ✓ El ecualizador tiene buena planicidad (ripple <= 3 dB)\n\n');
else
    fprintf('  ⚠ El ecualizador tiene ripple alto (> 3 dB)\n\n');
end

fprintf('=========================================================================\n');
fprintf('                  DISEÑO COMPLETADO EXITOSAMENTE\n');
fprintf('=========================================================================\n\n');
