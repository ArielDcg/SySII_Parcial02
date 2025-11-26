%% Diseño de Filtros Digitales IIR para 2.56kHz
% Segundo Parcial - Señales y Sistemas II
% Ecualizador de Audio de 4 bandas
%
% Este script diseña dos filtros digitales IIR pasabanda para la frecuencia
% de 2.56kHz usando:
%   1. Invarianza del impulso
%   2. Transformación bilineal
%
% Se comparan las respuestas en frecuencia de ambos métodos.

clear all;
close all;
clc;

%% Especificaciones del Sistema
fprintf('=========================================================================\n');
fprintf('         DISEÑO DE FILTROS DIGITALES IIR PARA 2.56 kHz\n');
fprintf('=========================================================================\n\n');

% Frecuencias centrales de las 4 bandas del ecualizador
f_central = [102.4, 512, 2560, 12800]; % Hz
fs = 64000; % Frecuencia de muestreo en Hz
T = 1/fs;   % Período de muestreo

% Banda asignada: 2.56kHz (tercera banda - PASABANDA)
fc = 2560; % Hz - Frecuencia central del filtro

fprintf('Frecuencia de muestreo: fs = %d Hz\n', fs);
fprintf('Período de muestreo: T = %.6f s\n', T);
fprintf('Frecuencia central del filtro: fc = %.0f Hz\n\n', fc);

%% Cálculo de Frecuencias Límite (igual que en el filtro analógico)
% Límite inferior (media geométrica entre 512Hz y 2560Hz)
f1_3dB = sqrt(512 * 2560);  % Frecuencia donde la atenuación es -3dB

% Límite superior (media geométrica entre 2560Hz y 12800Hz)
f2_3dB = sqrt(2560 * 12800);  % Frecuencia donde la atenuación es -3dB

% Frecuencias de rechazo (atenuación de -15dB)
f1_15dB = 512;    % Hz - Frecuencia central del filtro anterior
f2_15dB = 12800;  % Hz - Frecuencia central del filtro posterior

fprintf('Límite inferior de banda de paso (f1, -3dB): %.2f Hz\n', f1_3dB);
fprintf('Límite superior de banda de paso (f2, -3dB): %.2f Hz\n', f2_3dB);
fprintf('Frecuencia de rechazo inferior (f1s, -15dB): %.2f Hz\n', f1_15dB);
fprintf('Frecuencia de rechazo superior (f2s, -15dB): %.2f Hz\n\n', f2_15dB);

%% Conversión a Radianes/Segundo (Dominio Analógico)
w1 = 2*pi*f1_3dB;   % rad/s - Límite inferior de banda de paso
w2 = 2*pi*f2_3dB;   % rad/s - Límite superior de banda de paso
ws1 = 2*pi*f1_15dB; % rad/s - Frecuencia de rechazo inferior
ws2 = 2*pi*f2_15dB; % rad/s - Frecuencia de rechazo superior

% Frecuencia central del pasabanda (media geométrica en rad/s)
Wo = sqrt(w1*w2);  % rad/s
fc_calculada = Wo/(2*pi);

% Ancho de banda
BW = w2 - w1;  % rad/s

fprintf('Frecuencia central Wo: %.2f rad/s (%.2f Hz)\n', Wo, fc_calculada);
fprintf('Ancho de banda BW: %.2f rad/s (%.2f Hz)\n\n', BW, BW/(2*pi));

%% Especificaciones de Atenuación
Rp = 3;  % Atenuación máxima en la banda de paso (dB)
As = 15; % Atenuación mínima en la banda de rechazo (dB)

fprintf('Especificaciones:\n');
fprintf('  Rp (atenuación en banda de paso): %.2f dB\n', Rp);
fprintf('  As (atenuación en banda de rechazo): %.2f dB\n\n', As);

%% PASO 1: DISEÑO DEL FILTRO ANALÓGICO BUTTERWORTH
fprintf('=========================================================================\n');
fprintf('           PASO 1: DISEÑO DEL FILTRO ANALÓGICO DE REFERENCIA\n');
fprintf('=========================================================================\n\n');

% Cálculo del orden del filtro
[N_inferior, Wc_inferior] = buttord(w1, ws1, Rp, As, 's');
[N_superior, Wc_superior] = buttord(w2, ws2, Rp, As, 's');

% Seleccionar el orden mayor
N = max(N_inferior, N_superior);
fprintf('Orden del filtro prototipo pasabajos: N = %d\n', N);
fprintf('Orden del filtro pasabanda analógico: 2N = %d\n\n', 2*N);

% Diseño del filtro pasabajos prototipo
[b_lp, a_lp] = butter(N, 1, 's');

% Transformación a pasabanda
[b_bp_analog, a_bp_analog] = lp2bp(b_lp, a_lp, Wo, BW);

fprintf('Coeficientes del filtro pasabanda analógico:\n');
fprintf('Numerador (orden %d): ', length(b_bp_analog)-1);
fprintf('%.6f ', b_bp_analog);
fprintf('\n');
fprintf('Denominador (orden %d): ', length(a_bp_analog)-1);
fprintf('%.6f ', a_bp_analog);
fprintf('\n\n');

%% PASO 2: DISEÑO DEL FILTRO IIR CON INVARIANZA DEL IMPULSO
fprintf('=========================================================================\n');
fprintf('      PASO 2: DISEÑO DEL FILTRO IIR CON INVARIANZA DEL IMPULSO\n');
fprintf('=========================================================================\n\n');

% La función impinvar() convierte un filtro analógico a digital preservando
% la respuesta al impulso mediante:
%   H(z) = T * sum(residuos_i / (1 - exp(polos_i*T)*z^-1))
%
% IMPORTANTE: impinvar requiere que el filtro analógico tenga más polos que ceros
% para evitar aliasing en altas frecuencias.

fprintf('Aplicando invarianza del impulso...\n');
fprintf('Fórmula: h_d[n] = T * h_a(nT), donde T = 1/fs\n\n');

% Convertir el filtro analógico pasabanda a digital usando invarianza del impulso
[b_iir_imp, a_iir_imp] = impinvar(b_bp_analog, a_bp_analog, fs);

fprintf('Coeficientes del filtro IIR (Invarianza del Impulso):\n');
fprintf('Numerador (orden %d): ', length(b_iir_imp)-1);
fprintf('%.10f ', b_iir_imp);
fprintf('\n');
fprintf('Denominador (orden %d): ', length(a_iir_imp)-1);
fprintf('%.10f ', a_iir_imp);
fprintf('\n\n');

fprintf('Características del método de Invarianza del Impulso:\n');
fprintf('  + Preserva la forma de la respuesta al impulso del filtro analógico\n');
fprintf('  + Mapeo: s = (1/T) * ln(z)\n');
fprintf('  - Susceptible a aliasing en altas frecuencias\n');
fprintf('  - No preserva la respuesta en frecuencia en todo el rango\n\n');

%% PASO 3: DISEÑO DEL FILTRO IIR CON TRANSFORMACIÓN BILINEAL
fprintf('=========================================================================\n');
fprintf('       PASO 3: DISEÑO DEL FILTRO IIR CON TRANSFORMACIÓN BILINEAL\n');
fprintf('=========================================================================\n\n');

% La función bilinear() usa la transformación bilineal:
%   s = (2/T) * (z-1)/(z+1)  o equivalente:  s = 2*fs * (z-1)/(z+1)
%
% Esta transformación mapea todo el eje jω del plano s al círculo unitario
% en el plano z, evitando el aliasing pero introduciendo distorsión (warping)
% en frecuencia.

fprintf('Aplicando transformación bilineal...\n');
fprintf('Fórmula: s = (2/T) * (z-1)/(z+1) = %d * (z-1)/(z+1)\n\n', 2*fs);

% Convertir el filtro analógico pasabanda a digital usando transformación bilineal
[b_iir_bil, a_iir_bil] = bilinear(b_bp_analog, a_bp_analog, fs);

fprintf('Coeficientes del filtro IIR (Transformación Bilineal):\n');
fprintf('Numerador (orden %d): ', length(b_iir_bil)-1);
fprintf('%.10f ', b_iir_bil);
fprintf('\n');
fprintf('Denominador (orden %d): ', length(a_iir_bil)-1);
fprintf('%.10f ', a_iir_bil);
fprintf('\n\n');

fprintf('Características del método de Transformación Bilineal:\n');
fprintf('  + No presenta aliasing (todo el eje jω se mapea al círculo unitario)\n');
fprintf('  + Preserva la estabilidad del filtro analógico\n');
fprintf('  + Mapeo: s = (2/T) * (z-1)/(z+1)\n');
fprintf('  - Introduce distorsión en frecuencia (warping)\n');
fprintf('  - Relación no lineal: Ω = (2/T) * tan(ωT/2)\n\n');

%% PASO 4: CÁLCULO DE LAS RESPUESTAS EN FRECUENCIA
fprintf('=========================================================================\n');
fprintf('          PASO 4: CÁLCULO DE LAS RESPUESTAS EN FRECUENCIA\n');
fprintf('=========================================================================\n\n');

% Vector de frecuencias para análisis
N_points = 2048;
f_plot = linspace(0, fs/2, N_points); % De 0 a frecuencia de Nyquist
w_digital = 2*pi*f_plot/fs; % Frecuencia angular digital normalizada

% Respuesta en frecuencia del filtro analógico (para comparación)
w_analog = 2*pi*f_plot; % rad/s
[H_analog, w_a] = freqs(b_bp_analog, a_bp_analog, w_analog);
mag_analog_dB = 20*log10(abs(H_analog));
fase_analog_deg = angle(H_analog)*180/pi;

% Respuesta en frecuencia del filtro IIR con invarianza del impulso
[H_imp, w_imp] = freqz(b_iir_imp, a_iir_imp, w_digital);
mag_imp_dB = 20*log10(abs(H_imp));
fase_imp_deg = angle(H_imp)*180/pi;

% Respuesta en frecuencia del filtro IIR con transformación bilineal
[H_bil, w_bil] = freqz(b_iir_bil, a_iir_bil, w_digital);
mag_bil_dB = 20*log10(abs(H_bil));
fase_bil_deg = angle(H_bil)*180/pi;

fprintf('Respuestas en frecuencia calculadas para %d puntos.\n\n', N_points);

%% PASO 5: GRÁFICAS COMPARATIVAS
fprintf('=========================================================================\n');
fprintf('                   PASO 5: GRÁFICAS COMPARATIVAS\n');
fprintf('=========================================================================\n\n');

% Crear figura con gráficas comparativas
figure('Position', [50, 50, 1400, 900]);

% -------------------------------------------------------------------------
% Subplot 1: Comparación de Magnitud (escala logarítmica)
% -------------------------------------------------------------------------
subplot(3,2,1);
semilogx(f_plot, mag_analog_dB, 'k--', 'LineWidth', 2, 'DisplayName', 'Analógico (referencia)');
hold on;
semilogx(f_plot, mag_imp_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Invarianza Impulso');
semilogx(f_plot, mag_bil_dB, 'r-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Transf. Bilineal');
grid on;

% Marcar puntos importantes
plot(fc, 0, 'mo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'fc = 2560 Hz');
plot(f1_3dB, -3, 'go', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Límites -3dB');
plot(f2_3dB, -3, 'go', 'MarkerSize', 6, 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Líneas de referencia
yline(-3, '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
yline(-15, '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Comparación: Respuesta en Magnitud (Escala Log)');
legend('Location', 'southwest', 'FontSize', 8);
xlim([100 20000]);
ylim([-60 5]);

% -------------------------------------------------------------------------
% Subplot 2: Comparación de Magnitud (zoom en banda de paso)
% -------------------------------------------------------------------------
subplot(3,2,2);
plot(f_plot, mag_analog_dB, 'k--', 'LineWidth', 2, 'DisplayName', 'Analógico (referencia)');
hold on;
plot(f_plot, mag_imp_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Invarianza Impulso');
plot(f_plot, mag_bil_dB, 'r-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Transf. Bilineal');
grid on;

% Marcar puntos importantes
plot(fc, 0, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
plot(f1_3dB, -3, 'go', 'MarkerSize', 6, 'LineWidth', 1.5);
plot(f2_3dB, -3, 'go', 'MarkerSize', 6, 'LineWidth', 1.5);

yline(-3, '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Zoom en Banda de Paso');
legend('Location', 'southwest', 'FontSize', 8);
xlim([500 15000]);
ylim([-25 2]);

% -------------------------------------------------------------------------
% Subplot 3: Comparación de Fase (escala logarítmica)
% -------------------------------------------------------------------------
subplot(3,2,3);
semilogx(f_plot, fase_analog_deg, 'k--', 'LineWidth', 2, 'DisplayName', 'Analógico (referencia)');
hold on;
semilogx(f_plot, fase_imp_deg, 'b-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Invarianza Impulso');
semilogx(f_plot, fase_bil_deg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'IIR - Transf. Bilineal');
grid on;

% Marcar frecuencia central
xline(fc, '--m', 'LineWidth', 1);

xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
title('Comparación: Respuesta en Fase (Escala Log)');
legend('Location', 'southwest', 'FontSize', 8);
xlim([100 20000]);

% -------------------------------------------------------------------------
% Subplot 4: Diferencia de Magnitud respecto al filtro analógico
% -------------------------------------------------------------------------
subplot(3,2,4);
semilogx(f_plot, mag_imp_dB - mag_analog_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error Invarianza Impulso');
hold on;
semilogx(f_plot, mag_bil_dB - mag_analog_dB, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Error Transf. Bilineal');
grid on;
yline(0, '--k', 'LineWidth', 1);

xlabel('Frecuencia (Hz)');
ylabel('Diferencia (dB)');
title('Error de Magnitud respecto al Filtro Analógico');
legend('Location', 'best', 'FontSize', 8);
xlim([100 20000]);

% -------------------------------------------------------------------------
% Subplot 5: Diferencia de Fase respecto al filtro analógico
% -------------------------------------------------------------------------
subplot(3,2,5);
semilogx(f_plot, fase_imp_deg - fase_analog_deg, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Error Invarianza Impulso');
hold on;
semilogx(f_plot, fase_bil_deg - fase_analog_deg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Error Transf. Bilineal');
grid on;
yline(0, '--k', 'LineWidth', 1);

xlabel('Frecuencia (Hz)');
ylabel('Diferencia (grados)');
title('Error de Fase respecto al Filtro Analógico');
legend('Location', 'best', 'FontSize', 8);
xlim([100 20000]);

% -------------------------------------------------------------------------
% Subplot 6: Diagrama de polos y ceros
% -------------------------------------------------------------------------
subplot(3,2,6);
zplane(b_iir_imp, a_iir_imp);
title('Diagrama Polos-Ceros: Invarianza del Impulso');
grid on;

% Agregar título general
sgtitle('Comparación de Filtros IIR: Invarianza del Impulso vs Transformación Bilineal', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% PASO 6: GRÁFICA ADICIONAL - Diagrama de polos y ceros (Bilineal)
figure('Position', [100, 100, 800, 600]);
zplane(b_iir_bil, a_iir_bil);
title('Diagrama Polos-Ceros: Transformación Bilineal', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% VERIFICACIÓN DE ESPECIFICACIONES
fprintf('=========================================================================\n');
fprintf('              VERIFICACIÓN DE ESPECIFICACIONES\n');
fprintf('=========================================================================\n\n');

% Función auxiliar para encontrar el índice más cercano a una frecuencia
find_idx = @(f_vec, f_target) find(abs(f_vec - f_target) == min(abs(f_vec - f_target)), 1);

% Verificar invarianza del impulso
fprintf('--- FILTRO IIR CON INVARIANZA DEL IMPULSO ---\n');
idx_fc_imp = find_idx(f_plot, fc);
fprintf('Atenuación en fc = %.0f Hz: %.2f dB (esperado: ~0 dB)\n', fc, mag_imp_dB(idx_fc_imp));

idx_f1_imp = find_idx(f_plot, f1_3dB);
fprintf('Atenuación en f1 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f1_3dB, mag_imp_dB(idx_f1_imp));

idx_f2_imp = find_idx(f_plot, f2_3dB);
fprintf('Atenuación en f2 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f2_3dB, mag_imp_dB(idx_f2_imp));

idx_fs1_imp = find_idx(f_plot, f1_15dB);
fprintf('Atenuación en fs1 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f1_15dB, mag_imp_dB(idx_fs1_imp));

idx_fs2_imp = find_idx(f_plot, f2_15dB);
fprintf('Atenuación en fs2 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f2_15dB, mag_imp_dB(idx_fs2_imp));

% Verificar transformación bilineal
fprintf('\n--- FILTRO IIR CON TRANSFORMACIÓN BILINEAL ---\n');
idx_fc_bil = find_idx(f_plot, fc);
fprintf('Atenuación en fc = %.0f Hz: %.2f dB (esperado: ~0 dB)\n', fc, mag_bil_dB(idx_fc_bil));

idx_f1_bil = find_idx(f_plot, f1_3dB);
fprintf('Atenuación en f1 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f1_3dB, mag_bil_dB(idx_f1_bil));

idx_f2_bil = find_idx(f_plot, f2_3dB);
fprintf('Atenuación en f2 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f2_3dB, mag_bil_dB(idx_f2_bil));

idx_fs1_bil = find_idx(f_plot, f1_15dB);
fprintf('Atenuación en fs1 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f1_15dB, mag_bil_dB(idx_fs1_bil));

idx_fs2_bil = find_idx(f_plot, f2_15dB);
fprintf('Atenuación en fs2 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f2_15dB, mag_bil_dB(idx_fs2_bil));

%% ANÁLISIS Y COMPARACIÓN
fprintf('\n=========================================================================\n');
fprintf('                    ANÁLISIS Y COMPARACIÓN\n');
fprintf('=========================================================================\n\n');

fprintf('INVARIANZA DEL IMPULSO:\n');
fprintf('  • Preserva la respuesta al impulso del filtro analógico\n');
fprintf('  • Puede presentar aliasing en altas frecuencias\n');
fprintf('  • Mejor aproximación en bajas frecuencias\n');
fprintf('  • La respuesta en frecuencia puede diferir del filtro analógico\n');
fprintf('  • Mapeo: s = (1/T) * ln(z)\n\n');

fprintf('TRANSFORMACIÓN BILINEAL:\n');
fprintf('  • No presenta aliasing (mapeo completo al círculo unitario)\n');
fprintf('  • Introduce distorsión en frecuencia (warping)\n');
fprintf('  • La relación entre frecuencias analógicas y digitales es no lineal\n');
fprintf('  • Preserva mejor la respuesta en frecuencia en todo el rango\n');
fprintf('  • Mapeo: s = (2/T) * (z-1)/(z+1)\n\n');

fprintf('DIFERENCIAS OBSERVADAS:\n');
% Calcular error cuadrático medio en la banda de paso
idx_banda_paso = (f_plot >= f1_3dB) & (f_plot <= f2_3dB);
error_imp = sqrt(mean((mag_imp_dB(idx_banda_paso) - mag_analog_dB(idx_banda_paso)).^2));
error_bil = sqrt(mean((mag_bil_dB(idx_banda_paso) - mag_analog_dB(idx_banda_paso)).^2));

fprintf('  • Error RMS en banda de paso (Invarianza Impulso): %.4f dB\n', error_imp);
fprintf('  • Error RMS en banda de paso (Transf. Bilineal): %.4f dB\n', error_bil);

if error_bil < error_imp
    fprintf('\n>>> La Transformación Bilineal presenta MENOR error en la banda de paso.\n');
    filtro_recomendado = 'Transformación Bilineal';
else
    fprintf('\n>>> La Invarianza del Impulso presenta MENOR error en la banda de paso.\n');
    filtro_recomendado = 'Invarianza del Impulso';
end

%% GUARDAR RESULTADOS
fprintf('\n=========================================================================\n');
fprintf('                      GUARDAR RESULTADOS\n');
fprintf('=========================================================================\n\n');

% Guardar coeficientes de ambos filtros
save('filtros_iir_2560Hz.mat', ...
     'b_iir_imp', 'a_iir_imp', 'b_iir_bil', 'a_iir_bil', ...
     'b_bp_analog', 'a_bp_analog', 'N', 'fs', 'fc', ...
     'f1_3dB', 'f2_3dB', 'f1_15dB', 'f2_15dB', 'Rp', 'As', ...
     'filtro_recomendado');
fprintf('Coeficientes guardados en: filtros_iir_2560Hz.mat\n');

% Guardar figuras
saveas(figure(1), 'comparacion_filtros_iir_2560Hz.png');
fprintf('Gráfica comparativa guardada en: comparacion_filtros_iir_2560Hz.png\n');

saveas(figure(2), 'polos_ceros_bilineal_2560Hz.png');
fprintf('Diagrama polos-ceros guardado en: polos_ceros_bilineal_2560Hz.png\n');

%% RESUMEN FINAL
fprintf('\n=========================================================================\n');
fprintf('                          RESUMEN FINAL\n');
fprintf('=========================================================================\n');
fprintf('Filtro: Butterworth Pasabanda Digital IIR\n');
fprintf('Frecuencia de muestreo: fs = %d Hz\n', fs);
fprintf('Frecuencia central: fc = %.0f Hz\n', fc);
fprintf('Orden del filtro: %d\n', 2*N);
fprintf('Banda de paso (-3dB): %.2f Hz - %.2f Hz\n', f1_3dB, f2_3dB);
fprintf('Frecuencias de rechazo (-15dB): %.0f Hz y %.0f Hz\n', f1_15dB, f2_15dB);
fprintf('\nFiltro recomendado: %s\n', filtro_recomendado);
fprintf('=========================================================================\n');
