%% Diseño de Filtros FIR usando Ventanas para Ecualizador de Audio
% Segundo Parcial - Señales y Sistemas II
% Ecualizador de 4 bandas usando filtros FIR

clear all;
close all;
clc;

%% Especificaciones generales
fs = 64000; % Frecuencia de muestreo en Hz
Ts = 1/fs;  % Periodo de muestreo

% Frecuencias centrales de cada banda
f_central = [102.4, 512, 2560, 12800]; % Hz

% Atenuaciones especificadas
Ap = 6;  % Atenuación en banda de paso (dB) - especificación FIR
As = 15; % Atenuación en banda de rechazo (dB)

fprintf('=== DISEÑO DE FILTROS FIR USANDO VENTANAS ===\n\n');
fprintf('Frecuencia de muestreo: %.0f Hz\n', fs);
fprintf('Frecuencias centrales: %.1f Hz, %.0f Hz, %.2f kHz, %.1f kHz\n\n', ...
    f_central(1), f_central(2), f_central(3)/1000, f_central(4)/1000);

%% Cálculo de frecuencias de corte (media geométrica)
% Las frecuencias de corte están a -3dB y se calculan como la media
% geométrica entre frecuencias centrales adyacentes

fc = zeros(1, 3);
fc(1) = sqrt(f_central(1) * f_central(2));  % Entre 102.4 y 512 Hz
fc(2) = sqrt(f_central(2) * f_central(3));  % Entre 512 y 2560 Hz
fc(3) = sqrt(f_central(3) * f_central(4));  % Entre 2560 y 12800 Hz

fprintf('--- Frecuencias de corte (-3dB) ---\n');
fprintf('fc1 = %.2f Hz (entre %.1f y %.0f Hz)\n', fc(1), f_central(1), f_central(2));
fprintf('fc2 = %.2f Hz (entre %.0f y %.0f Hz)\n', fc(2), f_central(2), f_central(3));
fprintf('fc3 = %.2f Hz (entre %.0f y %.0f Hz)\n\n', fc(3), f_central(3), f_central(4));

%% Cálculo del ancho de banda de transición (Δf) para cada filtro

% Filtro 1: Pasabajos (102.4 Hz)
% Banda de paso: 0 - fc(1) Hz
% Banda de rechazo: desde f_central(2) Hz
delta_f1 = f_central(2) - fc(1);

% Filtro 2: Pasabanda (512 Hz)
% Banda de paso: fc(1) - fc(2) Hz
% Bandas de rechazo: hasta f_central(1) y desde f_central(3)
delta_f2_lower = fc(1) - f_central(1);
delta_f2_upper = f_central(3) - fc(2);
delta_f2 = min(delta_f2_lower, delta_f2_upper);

% Filtro 3: Pasabanda (2.56 kHz)
% Banda de paso: fc(2) - fc(3) Hz
% Bandas de rechazo: hasta f_central(2) y desde f_central(4)
delta_f3_lower = fc(2) - f_central(2);
delta_f3_upper = f_central(4) - fc(3);
delta_f3 = min(delta_f3_lower, delta_f3_upper);

% Filtro 4: Pasaaltos (12.8 kHz)
% Banda de paso: desde fc(3) Hz
% Banda de rechazo: hasta f_central(3)
delta_f4 = fc(3) - f_central(3);

fprintf('--- Anchos de banda de transición (Δf) ---\n');
fprintf('Filtro 1 (Pasabajos %.1f Hz):  Δf = %.2f Hz\n', f_central(1), delta_f1);
fprintf('Filtro 2 (Pasabanda %.0f Hz):  Δf = %.2f Hz\n', f_central(2), delta_f2);
fprintf('Filtro 3 (Pasabanda %.2f kHz): Δf = %.2f Hz\n', f_central(3)/1000, delta_f3);
fprintf('Filtro 4 (Pasaaltos %.1f kHz): Δf = %.2f Hz\n\n', f_central(4)/1000, delta_f4);

%% Selección de ventana
% Según las especificaciones: As = 15 dB
% Ventanas disponibles y sus atenuaciones:
% - Rectangular: 4π/M,  -13 dB (insuficiente)
% - Bartlett:    8π/M,  -25 dB (suficiente)
% - Hanning:     8π/M,  -31 dB (suficiente, mejor)
% - Hamming:     8π/M,  -41 dB (sobre-especificado)
% - Blackman:    12π/M, -57 dB (sobre-especificado)

% Usaremos Hanning por tener buen balance entre atenuación y ancho del lóbulo
ventana_tipo = 'hanning';
At_ventana = 31; % dB de atenuación de la ventana Hanning
ancho_lobulo = 8; % 8π/M para Hanning

fprintf('--- Selección de ventana ---\n');
fprintf('Ventana seleccionada: %s\n', upper(ventana_tipo));
fprintf('Atenuación del lóbulo secundario: %.0f dB\n', At_ventana);
fprintf('Ancho del lóbulo principal: %dπ/M\n\n', ancho_lobulo);

%% Cálculo del orden de los filtros usando la aproximación de Harris
% N = (fs * AT(dB)) / (Δf * 22)

N1 = ceil((fs * At_ventana) / (delta_f1 * 22));
N2 = ceil((fs * At_ventana) / (delta_f2 * 22));
N3 = ceil((fs * At_ventana) / (delta_f3 * 22));
N4 = ceil((fs * At_ventana) / (delta_f4 * 22));

% Asegurar que los órdenes sean impares para tener simetría
if mod(N1, 2) == 0, N1 = N1 + 1; end
if mod(N2, 2) == 0, N2 = N2 + 1; end
if mod(N3, 2) == 0, N3 = N3 + 1; end
if mod(N4, 2) == 0, N4 = N4 + 1; end

fprintf('--- Órdenes de los filtros (Fórmula de Harris) ---\n');
fprintf('Filtro 1 (Pasabajos):  N = %d (M = %d coeficientes)\n', N1, N1+1);
fprintf('Filtro 2 (Pasabanda):  N = %d (M = %d coeficientes)\n', N2, N2+1);
fprintf('Filtro 3 (Pasabanda):  N = %d (M = %d coeficientes)\n', N3, N3+1);
fprintf('Filtro 4 (Pasaaltos):  N = %d (M = %d coeficientes)\n\n', N4, N4+1);

%% Diseño de los filtros FIR usando el método de ventanas

fprintf('--- Diseñando filtros FIR ---\n');

% Normalizar frecuencias de corte (respecto a frecuencia de Nyquist)
wc1 = fc(1) / (fs/2);
wc2 = [fc(1), fc(2)] / (fs/2);
wc3 = [fc(2), fc(3)] / (fs/2);
wc4 = fc(3) / (fs/2);

% Filtro 1: Pasabajos
fprintf('Diseñando Filtro 1 (Pasabajos %.1f Hz)...\n', f_central(1));
b1 = fir1(N1, wc1, 'low', hanning(N1+1));
fprintf('  Coeficientes: %d\n', length(b1));

% Filtro 2: Pasabanda
fprintf('Diseñando Filtro 2 (Pasabanda %.0f Hz)...\n', f_central(2));
b2 = fir1(N2, wc2, 'bandpass', hanning(N2+1));
fprintf('  Coeficientes: %d\n', length(b2));

% Filtro 3: Pasabanda
fprintf('Diseñando Filtro 3 (Pasabanda %.2f kHz)...\n', f_central(3)/1000);
b3 = fir1(N3, wc3, 'bandpass', hanning(N3+1));
fprintf('  Coeficientes: %d\n', length(b3));

% Filtro 4: Pasaaltos
fprintf('Diseñando Filtro 4 (Pasaaltos %.1f kHz)...\n', f_central(4)/1000);
b4 = fir1(N4, wc4, 'high', hanning(N4+1));
fprintf('  Coeficientes: %d\n\n', length(b4));

%% Guardar los filtros diseñados
save('filtros_FIR_diseñados.mat', 'b1', 'b2', 'b3', 'b4', 'fs', ...
     'f_central', 'fc', 'N1', 'N2', 'N3', 'N4', 'ventana_tipo');

fprintf('Filtros guardados en: filtros_FIR_diseñados.mat\n\n');

%% Calcular y graficar respuestas en frecuencia
fprintf('Generando gráficas de respuesta en frecuencia...\n');

% Número de puntos para el análisis de frecuencia
Nfft = 8192;
f_axis = linspace(0, fs/2, Nfft/2+1);

% Respuestas en frecuencia
[H1, w1] = freqz(b1, 1, Nfft, fs);
[H2, w2] = freqz(b2, 1, Nfft, fs);
[H3, w3] = freqz(b3, 1, Nfft, fs);
[H4, w4] = freqz(b4, 1, Nfft, fs);

% Magnitud en dB
H1_dB = 20*log10(abs(H1));
H2_dB = 20*log10(abs(H2));
H3_dB = 20*log10(abs(H3));
H4_dB = 20*log10(abs(H4));

% Fase en grados
H1_phase = angle(H1) * 180/pi;
H2_phase = angle(H2) * 180/pi;
H3_phase = angle(H3) * 180/pi;
H4_phase = angle(H4) * 180/pi;

%% Figura 1: Respuestas en Magnitud
figure('Name', 'Respuestas en Magnitud - Filtros FIR', 'Position', [100 100 1200 800]);

% Filtro 1
subplot(4,1,1);
semilogx(w1, H1_dB, 'b', 'LineWidth', 1.5);
grid on;
hold on;
% Marcar frecuencia central y de corte
xline(f_central(1), 'r--', sprintf('f_c=%.1f Hz', f_central(1)), 'LineWidth', 1);
xline(fc(1), 'g--', sprintf('f_{-3dB}=%.1f Hz', fc(1)), 'LineWidth', 1);
xline(f_central(2), 'm--', sprintf('f_{-15dB}=%.0f Hz', f_central(2)), 'LineWidth', 1);
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Filtro 1: Pasabajos (%.1f Hz) - Orden N=%d', f_central(1), N1));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([10 fs/2]);
legend('H_1(f)', 'Location', 'southwest');

% Filtro 2
subplot(4,1,2);
semilogx(w2, H2_dB, 'b', 'LineWidth', 1.5);
grid on;
hold on;
xline(f_central(2), 'r--', sprintf('f_c=%.0f Hz', f_central(2)), 'LineWidth', 1);
xline(fc(1), 'g--', sprintf('f_{-3dB,L}=%.1f Hz', fc(1)), 'LineWidth', 1);
xline(fc(2), 'g--', sprintf('f_{-3dB,H}=%.1f Hz', fc(2)), 'LineWidth', 1);
xline(f_central(1), 'm--', sprintf('f_{-15dB,L}=%.1f Hz', f_central(1)), 'LineWidth', 1);
xline(f_central(3), 'm--', sprintf('f_{-15dB,H}=%.0f Hz', f_central(3)), 'LineWidth', 1);
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Filtro 2: Pasabanda (%.0f Hz) - Orden N=%d', f_central(2), N2));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([10 fs/2]);
legend('H_2(f)', 'Location', 'southwest');

% Filtro 3
subplot(4,1,3);
semilogx(w3, H3_dB, 'b', 'LineWidth', 1.5);
grid on;
hold on;
xline(f_central(3), 'r--', sprintf('f_c=%.2f kHz', f_central(3)/1000), 'LineWidth', 1);
xline(fc(2), 'g--', sprintf('f_{-3dB,L}=%.1f Hz', fc(2)), 'LineWidth', 1);
xline(fc(3), 'g--', sprintf('f_{-3dB,H}=%.1f kHz', fc(3)/1000), 'LineWidth', 1);
xline(f_central(2), 'm--', sprintf('f_{-15dB,L}=%.0f Hz', f_central(2)), 'LineWidth', 1);
xline(f_central(4), 'm--', sprintf('f_{-15dB,H}=%.1f kHz', f_central(4)/1000), 'LineWidth', 1);
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Filtro 3: Pasabanda (%.2f kHz) - Orden N=%d', f_central(3)/1000, N3));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([10 fs/2]);
legend('H_3(f)', 'Location', 'southwest');

% Filtro 4
subplot(4,1,4);
semilogx(w4, H4_dB, 'b', 'LineWidth', 1.5);
grid on;
hold on;
xline(f_central(4), 'r--', sprintf('f_c=%.1f kHz', f_central(4)/1000), 'LineWidth', 1);
xline(fc(3), 'g--', sprintf('f_{-3dB}=%.1f kHz', fc(3)/1000), 'LineWidth', 1);
xline(f_central(3), 'm--', sprintf('f_{-15dB}=%.2f kHz', f_central(3)/1000), 'LineWidth', 1);
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Filtro 4: Pasaaltos (%.1f kHz) - Orden N=%d', f_central(4)/1000, N4));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([10 fs/2]);
legend('H_4(f)', 'Location', 'southwest');

saveas(gcf, 'FIR_Respuestas_Magnitud.png');

%% Figura 2: Respuestas en Fase
figure('Name', 'Respuestas en Fase - Filtros FIR', 'Position', [150 150 1200 800]);

% Filtro 1
subplot(4,1,1);
semilogx(w1, H1_phase, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 1: Pasabajos (%.1f Hz) - Fase', f_central(1)));
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([10 fs/2]);

% Filtro 2
subplot(4,1,2);
semilogx(w2, H2_phase, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 2: Pasabanda (%.0f Hz) - Fase', f_central(2)));
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([10 fs/2]);

% Filtro 3
subplot(4,1,3);
semilogx(w3, H3_phase, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 3: Pasabanda (%.2f kHz) - Fase', f_central(3)/1000));
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([10 fs/2]);

% Filtro 4
subplot(4,1,4);
semilogx(w4, H4_phase, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 4: Pasaaltos (%.1f kHz) - Fase', f_central(4)/1000));
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([10 fs/2]);

saveas(gcf, 'FIR_Respuestas_Fase.png');

%% Figura 3: Respuesta combinada de todos los filtros
figure('Name', 'Respuestas Combinadas - Todos los Filtros FIR', 'Position', [200 200 1000 600]);

subplot(2,1,1);
semilogx(w1, H1_dB, 'b', 'LineWidth', 1.5); hold on;
semilogx(w2, H2_dB, 'r', 'LineWidth', 1.5);
semilogx(w3, H3_dB, 'g', 'LineWidth', 1.5);
semilogx(w4, H4_dB, 'm', 'LineWidth', 1.5);
grid on;
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title('Respuestas en Magnitud - Todos los Filtros FIR');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([10 fs/2]);
legend('Pasabajos (102.4 Hz)', 'Pasabanda (512 Hz)', ...
       'Pasabanda (2.56 kHz)', 'Pasaaltos (12.8 kHz)', 'Location', 'southwest');

subplot(2,1,2);
semilogx(w1, H1_phase, 'b', 'LineWidth', 1.5); hold on;
semilogx(w2, H2_phase, 'r', 'LineWidth', 1.5);
semilogx(w3, H3_phase, 'g', 'LineWidth', 1.5);
semilogx(w4, H4_phase, 'm', 'LineWidth', 1.5);
grid on;
title('Respuestas en Fase - Todos los Filtros FIR');
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([10 fs/2]);
legend('Pasabajos (102.4 Hz)', 'Pasabanda (512 Hz)', ...
       'Pasabanda (2.56 kHz)', 'Pasaaltos (12.8 kHz)', 'Location', 'southwest');

saveas(gcf, 'FIR_Respuestas_Combinadas.png');

%% Verificación de especificaciones
fprintf('--- Verificación de especificaciones ---\n\n');

% Función auxiliar para encontrar ganancia a una frecuencia específica
find_gain = @(H, w, f_target) interp1(w, 20*log10(abs(H)), f_target);

% Filtro 1 - Pasabajos
fprintf('Filtro 1 (Pasabajos %.1f Hz):\n', f_central(1));
fprintf('  Ganancia a %.1f Hz (central):    %.2f dB\n', f_central(1), find_gain(H1, w1, f_central(1)));
fprintf('  Ganancia a %.2f Hz (corte):      %.2f dB\n', fc(1), find_gain(H1, w1, fc(1)));
fprintf('  Ganancia a %.0f Hz (rechazo):    %.2f dB\n\n', f_central(2), find_gain(H1, w1, f_central(2)));

% Filtro 2 - Pasabanda
fprintf('Filtro 2 (Pasabanda %.0f Hz):\n', f_central(2));
fprintf('  Ganancia a %.0f Hz (central):    %.2f dB\n', f_central(2), find_gain(H2, w2, f_central(2)));
fprintf('  Ganancia a %.2f Hz (corte inf):  %.2f dB\n', fc(1), find_gain(H2, w2, fc(1)));
fprintf('  Ganancia a %.2f Hz (corte sup):  %.2f dB\n', fc(2), find_gain(H2, w2, fc(2)));
fprintf('  Ganancia a %.1f Hz (rechazo):    %.2f dB\n', f_central(1), find_gain(H2, w2, f_central(1)));
fprintf('  Ganancia a %.0f Hz (rechazo):    %.2f dB\n\n', f_central(3), find_gain(H2, w2, f_central(3)));

% Filtro 3 - Pasabanda
fprintf('Filtro 3 (Pasabanda %.2f kHz):\n', f_central(3)/1000);
fprintf('  Ganancia a %.0f Hz (central):    %.2f dB\n', f_central(3), find_gain(H3, w3, f_central(3)));
fprintf('  Ganancia a %.2f Hz (corte inf):  %.2f dB\n', fc(2), find_gain(H3, w3, fc(2)));
fprintf('  Ganancia a %.2f Hz (corte sup):  %.2f dB\n', fc(3), find_gain(H3, w3, fc(3)));
fprintf('  Ganancia a %.0f Hz (rechazo):    %.2f dB\n', f_central(2), find_gain(H3, w3, f_central(2)));
fprintf('  Ganancia a %.0f Hz (rechazo):    %.2f dB\n\n', f_central(4), find_gain(H3, w3, f_central(4)));

% Filtro 4 - Pasaaltos
fprintf('Filtro 4 (Pasaaltos %.1f kHz):\n', f_central(4)/1000);
fprintf('  Ganancia a %.0f Hz (central):    %.2f dB\n', f_central(4), find_gain(H4, w4, f_central(4)));
fprintf('  Ganancia a %.2f Hz (corte):      %.2f dB\n', fc(3), find_gain(H4, w4, fc(3)));
fprintf('  Ganancia a %.0f Hz (rechazo):    %.2f dB\n\n', f_central(3), find_gain(H4, w4, f_central(3)));

fprintf('=== DISEÑO COMPLETADO ===\n');
fprintf('Los filtros FIR han sido diseñados exitosamente.\n');
fprintf('Archivos generados:\n');
fprintf('  - filtros_FIR_diseñados.mat\n');
fprintf('  - FIR_Respuestas_Magnitud.png\n');
fprintf('  - FIR_Respuestas_Fase.png\n');
fprintf('  - FIR_Respuestas_Combinadas.png\n');
