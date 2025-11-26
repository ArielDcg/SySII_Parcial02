%% Diseño de Filtro Pasabanda FIR 2.56 kHz usando Ventanas
% Segundo Parcial - Señales y Sistemas II
% Comparación entre ventanas Hanning y Hamming

clear all;
close all;
clc;

%% Especificaciones del sistema
fs = 64000; % Frecuencia de muestreo en Hz

% Frecuencias centrales del ecualizador
f_central = [102.4, 512, 2560, 12800]; % Hz

% Atenuaciones especificadas
Ap = 6;  % Atenuación en banda de paso (dB) - especificación FIR
As = 15; % Atenuación en banda de rechazo (dB)

fprintf('=== DISEÑO DE FILTRO PASABANDA FIR 2.56 kHz ===\n');
fprintf('Comparación de ventanas: HANNING vs HAMMING\n\n');

%% Especificaciones del filtro pasabanda 2.56 kHz
% Frecuencia central: 2560 Hz
fc_target = 2560; % Hz

% Cálculo de frecuencias de corte (media geométrica)
fc_lower = sqrt(f_central(2) * f_central(3)); % Entre 512 Hz y 2560 Hz
fc_upper = sqrt(f_central(3) * f_central(4)); % Entre 2560 Hz y 12800 Hz

fprintf('--- Especificaciones del filtro ---\n');
fprintf('Frecuencia central:        %.0f Hz\n', fc_target);
fprintf('Frecuencia de corte inferior (-3dB): %.2f Hz\n', fc_lower);
fprintf('Frecuencia de corte superior (-3dB): %.2f Hz\n', fc_upper);
fprintf('Frecuencia de muestreo:    %.0f Hz\n', fs);
fprintf('Atenuación en banda de paso:  %.0f dB\n', Ap);
fprintf('Atenuación en banda de rechazo: %.0f dB\n\n', As);

%% Cálculo del ancho de banda de transición
% Para pasabanda:
% Banda de paso: fc_lower - fc_upper Hz (-3dB)
% Bandas de rechazo: hasta f_central(2) y desde f_central(4) (-15dB)

delta_f_lower = fc_lower - f_central(2);
delta_f_upper = f_central(4) - fc_upper;
delta_f = min(delta_f_lower, delta_f_upper);

fprintf('--- Ancho de banda de transición ---\n');
fprintf('Δf inferior: %.2f Hz (desde %.0f Hz hasta %.2f Hz)\n', ...
        delta_f_lower, f_central(2), fc_lower);
fprintf('Δf superior: %.2f Hz (desde %.2f Hz hasta %.0f Hz)\n', ...
        delta_f_upper, fc_upper, f_central(4));
fprintf('Δf crítico (mínimo): %.2f Hz\n\n', delta_f);

%% Propiedades de las ventanas
fprintf('--- Propiedades de las ventanas ---\n\n');

% VENTANA HANNING
fprintf('VENTANA HANNING:\n');
fprintf('  Ancho del lóbulo principal: 8π/M\n');
fprintf('  Pico del lóbulo secundario: -31 dB\n');
At_hanning = 31; % dB
ancho_hanning = 8; % 8π/M

% Cálculo del orden usando fórmula de Harris
% N = (fs * AT(dB)) / (Δf * 22)
N_hanning = ceil((fs * At_hanning) / (delta_f * 22));
if mod(N_hanning, 2) == 0
    N_hanning = N_hanning + 1; % Hacer impar para simetría
end
M_hanning = N_hanning + 1; % Número de coeficientes
fprintf('  Orden calculado (Harris): N = %d\n', N_hanning);
fprintf('  Número de coeficientes: M = %d\n\n', M_hanning);

% VENTANA HAMMING
fprintf('VENTANA HAMMING:\n');
fprintf('  Ancho del lóbulo principal: 8π/M\n');
fprintf('  Pico del lóbulo secundario: -41 dB\n');
At_hamming = 41; % dB
ancho_hamming = 8; % 8π/M

% Cálculo del orden
N_hamming = ceil((fs * At_hamming) / (delta_f * 22));
if mod(N_hamming, 2) == 0
    N_hamming = N_hamming + 1;
end
M_hamming = N_hamming + 1;
fprintf('  Orden calculado (Harris): N = %d\n', N_hamming);
fprintf('  Número de coeficientes: M = %d\n\n', M_hamming);

%% Diseño de los filtros usando ambas ventanas
fprintf('--- Diseñando filtros ---\n');

% Normalizar frecuencias de corte
wc = [fc_lower, fc_upper] / (fs/2);

% Filtro con ventana Hanning
fprintf('Diseñando filtro con ventana Hanning...\n');
b_hanning = fir1(N_hanning, wc, 'bandpass', hanning(M_hanning));
fprintf('  Filtro diseñado: %d coeficientes\n', length(b_hanning));

% Filtro con ventana Hamming
fprintf('Diseñando filtro con ventana Hamming...\n');
b_hamming = fir1(N_hamming, wc, 'bandpass', hamming(M_hamming));
fprintf('  Filtro diseñado: %d coeficientes\n\n', length(b_hamming));

%% Análisis de respuesta en frecuencia
fprintf('--- Análisis de respuesta en frecuencia ---\n');

Nfft = 8192;
[H_hanning, w_hanning] = freqz(b_hanning, 1, Nfft, fs);
[H_hamming, w_hamming] = freqz(b_hamming, 1, Nfft, fs);

% Magnitud en dB
H_hanning_dB = 20*log10(abs(H_hanning));
H_hamming_dB = 20*log10(abs(H_hamming));

% Fase en grados
H_hanning_phase = angle(H_hanning) * 180/pi;
H_hamming_phase = angle(H_hamming) * 180/pi;

% Retardo de grupo
[gd_hanning, w_gd_hanning] = grpdelay(b_hanning, 1, Nfft, fs);
[gd_hamming, w_gd_hamming] = grpdelay(b_hamming, 1, Nfft, fs);

%% Verificación de especificaciones
fprintf('\nVERIFICACIÓN DE ESPECIFICACIONES:\n\n');

% Función auxiliar para encontrar ganancia a una frecuencia específica
find_gain = @(H, w, f_target) interp1(w, 20*log10(abs(H)), f_target);

% VENTANA HANNING
fprintf('FILTRO CON VENTANA HANNING:\n');
gain_center_han = find_gain(H_hanning, w_hanning, fc_target);
gain_lower_han = find_gain(H_hanning, w_hanning, fc_lower);
gain_upper_han = find_gain(H_hanning, w_hanning, fc_upper);
gain_rej_low_han = find_gain(H_hanning, w_hanning, f_central(2));
gain_rej_high_han = find_gain(H_hanning, w_hanning, f_central(4));

fprintf('  Ganancia a %.0f Hz (central):      %6.2f dB (esperado: ~0 dB)\n', ...
        fc_target, gain_center_han);
fprintf('  Ganancia a %.2f Hz (corte inf):    %6.2f dB (esperado: -3 dB)\n', ...
        fc_lower, gain_lower_han);
fprintf('  Ganancia a %.2f Hz (corte sup):    %6.2f dB (esperado: -3 dB)\n', ...
        fc_upper, gain_upper_han);
fprintf('  Ganancia a %.0f Hz (rechazo inf):  %6.2f dB (esperado: <-15 dB)\n', ...
        f_central(2), gain_rej_low_han);
fprintf('  Ganancia a %.0f Hz (rechazo sup):  %6.2f dB (esperado: <-15 dB)\n', ...
        f_central(4), gain_rej_high_han);

% Verificar cumplimiento
cumple_han = (abs(gain_center_han) <= Ap) && ...
             (abs(gain_lower_han + 3) <= 1) && ...
             (abs(gain_upper_han + 3) <= 1) && ...
             (gain_rej_low_han <= -As) && ...
             (gain_rej_high_han <= -As);

if cumple_han
    fprintf('  ✓ CUMPLE con las especificaciones\n\n');
else
    fprintf('  ✗ NO CUMPLE con todas las especificaciones\n\n');
end

% VENTANA HAMMING
fprintf('FILTRO CON VENTANA HAMMING:\n');
gain_center_ham = find_gain(H_hamming, w_hamming, fc_target);
gain_lower_ham = find_gain(H_hamming, w_hamming, fc_lower);
gain_upper_ham = find_gain(H_hamming, w_hamming, fc_upper);
gain_rej_low_ham = find_gain(H_hamming, w_hamming, f_central(2));
gain_rej_high_ham = find_gain(H_hamming, w_hamming, f_central(4));

fprintf('  Ganancia a %.0f Hz (central):      %6.2f dB (esperado: ~0 dB)\n', ...
        fc_target, gain_center_ham);
fprintf('  Ganancia a %.2f Hz (corte inf):    %6.2f dB (esperado: -3 dB)\n', ...
        fc_lower, gain_lower_ham);
fprintf('  Ganancia a %.2f Hz (corte sup):    %6.2f dB (esperado: -3 dB)\n', ...
        fc_upper, gain_upper_ham);
fprintf('  Ganancia a %.0f Hz (rechazo inf):  %6.2f dB (esperado: <-15 dB)\n', ...
        f_central(2), gain_rej_low_ham);
fprintf('  Ganancia a %.0f Hz (rechazo sup):  %6.2f dB (esperado: <-15 dB)\n', ...
        f_central(4), gain_rej_high_ham);

cumple_ham = (abs(gain_center_ham) <= Ap) && ...
             (abs(gain_lower_ham + 3) <= 1) && ...
             (abs(gain_upper_ham + 3) <= 1) && ...
             (gain_rej_low_ham <= -As) && ...
             (gain_rej_high_ham <= -As);

if cumple_ham
    fprintf('  ✓ CUMPLE con las especificaciones\n\n');
else
    fprintf('  ✗ NO CUMPLE con todas las especificaciones\n\n');
end

%% Comparación entre ventanas
fprintf('--- COMPARACIÓN ENTRE VENTANAS ---\n\n');
fprintf('Parámetro                  | Hanning      | Hamming\n');
fprintf('--------------------------------------------------------\n');
fprintf('Orden del filtro (N)       | %-12d | %d\n', N_hanning, N_hamming);
fprintf('Núm. de coeficientes (M)   | %-12d | %d\n', M_hanning, M_hamming);
fprintf('Atenuación lóbulo sec.     | %-12.0f | %.0f dB\n', At_hanning, At_hamming);
fprintf('Ganancia en frecuencia central | %.2f         | %.2f dB\n', gain_center_han, gain_center_ham);
fprintf('Ganancia en fc inferior    | %.2f         | %.2f dB\n', gain_lower_han, gain_lower_ham);
fprintf('Ganancia en fc superior    | %.2f         | %.2f dB\n', gain_upper_han, gain_upper_ham);
fprintf('Atenuación en 512 Hz       | %.2f        | %.2f dB\n', gain_rej_low_han, gain_rej_low_ham);
fprintf('Atenuación en 12.8 kHz     | %.2f        | %.2f dB\n', gain_rej_high_han, gain_rej_high_ham);
fprintf('Retardo (muestras)         | %-12.1f | %.1f\n', N_hanning/2, N_hamming/2);
fprintf('Retardo (ms)               | %-12.3f | %.3f\n', N_hanning/(2*fs)*1000, N_hamming/(2*fs)*1000);
fprintf('\n');

%% Gráficas comparativas
% Figura 1: Respuestas en Magnitud
figure('Name', 'Comparación Magnitud - Hanning vs Hamming', ...
       'Position', [100 100 1200 800]);

subplot(2,2,1);
semilogx(w_hanning, H_hanning_dB, 'b', 'LineWidth', 2);
grid on;
hold on;
xline(fc_target, 'r--', sprintf('f_c=%.0f Hz', fc_target), 'LineWidth', 1);
xline(fc_lower, 'g--', sprintf('f_{-3dB,L}=%.0f Hz', fc_lower), 'LineWidth', 1);
xline(fc_upper, 'g--', sprintf('f_{-3dB,H}=%.0f Hz', fc_upper), 'LineWidth', 1);
xline(f_central(2), 'm--', sprintf('%.0f Hz', f_central(2)), 'LineWidth', 1);
xline(f_central(4), 'm--', sprintf('%.0f Hz', f_central(4)), 'LineWidth', 1);
yline(0, 'k:', '0dB');
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Magnitud - Ventana Hanning (N=%d)', N_hanning));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([50 fs/2]);

subplot(2,2,2);
semilogx(w_hamming, H_hamming_dB, 'r', 'LineWidth', 2);
grid on;
hold on;
xline(fc_target, 'r--', sprintf('f_c=%.0f Hz', fc_target), 'LineWidth', 1);
xline(fc_lower, 'g--', sprintf('f_{-3dB,L}=%.0f Hz', fc_lower), 'LineWidth', 1);
xline(fc_upper, 'g--', sprintf('f_{-3dB,H}=%.0f Hz', fc_upper), 'LineWidth', 1);
xline(f_central(2), 'm--', sprintf('%.0f Hz', f_central(2)), 'LineWidth', 1);
xline(f_central(4), 'm--', sprintf('%.0f Hz', f_central(4)), 'LineWidth', 1);
yline(0, 'k:', '0dB');
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title(sprintf('Magnitud - Ventana Hamming (N=%d)', N_hamming));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([50 fs/2]);

subplot(2,2,3);
semilogx(w_hanning, H_hanning_dB, 'b', 'LineWidth', 2);
hold on;
semilogx(w_hamming, H_hamming_dB, 'r', 'LineWidth', 2);
grid on;
xline(fc_target, 'k--', 'LineWidth', 0.5);
xline(fc_lower, 'k--', 'LineWidth', 0.5);
xline(fc_upper, 'k--', 'LineWidth', 0.5);
yline(0, 'k:', '0dB');
yline(-3, 'k:', '-3dB');
yline(-15, 'k:', '-15dB');
title('Comparación de Magnitud');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
ylim([-80 5]);
xlim([50 fs/2]);
legend('Hanning', 'Hamming', 'Location', 'southwest');

subplot(2,2,4);
plot(w_hanning, H_hanning_dB, 'b', 'LineWidth', 2);
hold on;
plot(w_hamming, H_hamming_dB, 'r', 'LineWidth', 2);
grid on;
xline(fc_lower, 'k--', 'LineWidth', 0.5);
xline(fc_upper, 'k--', 'LineWidth', 0.5);
yline(0, 'k:', '0dB');
yline(-3, 'k:', '-3dB');
title('Detalle en Banda de Paso');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([fc_lower*0.5 fc_upper*1.5]);
ylim([-10 2]);
legend('Hanning', 'Hamming', 'Location', 'southwest');

saveas(gcf, 'FIR_2560Hz_Comparacion_Magnitud.png');

% Figura 2: Respuestas en Fase
figure('Name', 'Comparación Fase - Hanning vs Hamming', ...
       'Position', [150 150 1200 600]);

subplot(2,1,1);
semilogx(w_hanning, H_hanning_phase, 'b', 'LineWidth', 2);
hold on;
semilogx(w_hamming, H_hamming_phase, 'r', 'LineWidth', 2);
grid on;
title('Respuesta en Fase');
xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
xlim([50 fs/2]);
legend('Hanning', 'Hamming', 'Location', 'southwest');

subplot(2,1,2);
semilogx(w_gd_hanning, gd_hanning, 'b', 'LineWidth', 2);
hold on;
semilogx(w_gd_hamming, gd_hamming, 'r', 'LineWidth', 2);
grid on;
title('Retardo de Grupo');
xlabel('Frecuencia (Hz)');
ylabel('Retardo (muestras)');
xlim([50 fs/2]);
legend('Hanning', 'Hamming', 'Location', 'northeast');

saveas(gcf, 'FIR_2560Hz_Comparacion_Fase.png');

% Figura 3: Respuesta al impulso
figure('Name', 'Respuesta al Impulso - Hanning vs Hamming', ...
       'Position', [200 200 1200 600]);

subplot(2,1,1);
stem(0:length(b_hanning)-1, b_hanning, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Respuesta al Impulso - Ventana Hanning (M=%d coeficientes)', M_hanning));
xlabel('n (muestras)');
ylabel('h[n]');

subplot(2,1,2);
stem(0:length(b_hamming)-1, b_hamming, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Respuesta al Impulso - Ventana Hamming (M=%d coeficientes)', M_hamming));
xlabel('n (muestras)');
ylabel('h[n]');

saveas(gcf, 'FIR_2560Hz_Respuesta_Impulso.png');

% Figura 4: Forma de las ventanas
figure('Name', 'Ventanas Utilizadas', 'Position', [250 250 1000 600]);

n_han = 0:N_hanning;
w_han_func = hanning(M_hanning);

n_ham = 0:N_hamming;
w_ham_func = hamming(M_hamming);

subplot(2,1,1);
stem(n_han, w_han_func, 'b', 'LineWidth', 1.5);
grid on;
title(sprintf('Ventana Hanning (M=%d)', M_hanning));
xlabel('n (muestras)');
ylabel('w[n]');
ylim([0 1.1]);

subplot(2,1,2);
stem(n_ham, w_ham_func, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Ventana Hamming (M=%d)', M_hamming));
xlabel('n (muestras)');
ylabel('w[n]');
ylim([0 1.1]);

saveas(gcf, 'FIR_2560Hz_Ventanas.png');

%% Guardar resultados
save('filtro_FIR_2560Hz_comparacion.mat', ...
     'b_hanning', 'b_hamming', 'N_hanning', 'N_hamming', ...
     'H_hanning', 'H_hamming', 'w_hanning', 'w_hamming', ...
     'fs', 'fc_target', 'fc_lower', 'fc_upper', ...
     'gain_center_han', 'gain_center_ham', ...
     'gain_rej_low_han', 'gain_rej_low_ham', ...
     'gain_rej_high_han', 'gain_rej_high_ham');

fprintf('=== DISEÑO COMPLETADO ===\n');
fprintf('Archivos generados:\n');
fprintf('  - filtro_FIR_2560Hz_comparacion.mat\n');
fprintf('  - FIR_2560Hz_Comparacion_Magnitud.png\n');
fprintf('  - FIR_2560Hz_Comparacion_Fase.png\n');
fprintf('  - FIR_2560Hz_Respuesta_Impulso.png\n');
fprintf('  - FIR_2560Hz_Ventanas.png\n\n');

%% Recomendación
fprintf('--- RECOMENDACIÓN ---\n');
if N_hanning < N_hamming
    fprintf('La ventana HANNING requiere un orden menor (%d vs %d).\n', N_hanning, N_hamming);
else
    fprintf('La ventana HAMMING requiere un orden menor (%d vs %d).\n', N_hamming, N_hanning);
end

if cumple_han && cumple_ham
    fprintf('Ambas ventanas cumplen con las especificaciones.\n');
    fprintf('HANNING es más eficiente computacionalmente (menor orden).\n');
    fprintf('HAMMING ofrece mayor atenuación en banda de rechazo.\n\n');
    fprintf('Recomendación: HANNING para eficiencia, HAMMING para mayor rechazo.\n');
elseif cumple_han
    fprintf('Solo la ventana HANNING cumple con las especificaciones.\n');
    fprintf('Recomendación: Usar HANNING.\n');
elseif cumple_ham
    fprintf('Solo la ventana HAMMING cumple con las especificaciones.\n');
    fprintf('Recomendación: Usar HAMMING.\n');
else
    fprintf('Ninguna ventana cumple completamente. Se requiere ajuste.\n');
end
