%% Diseño de Filtro Butterworth Pasabanda para 2.56kHz
% Segundo Parcial - Señales y Sistemas II
% Ecualizador de Audio de 4 bandas
%
% Este script diseña un filtro Butterworth pasabanda analógico para la
% frecuencia de 2.56kHz usando las especificaciones del parcial.

clear all;
close all;
clc;

%% Especificaciones del Sistema
% Frecuencias centrales de las 4 bandas del ecualizador
f_central = [102.4, 512, 2560, 12800]; % Hz
fs = 64000; % Frecuencia de muestreo en Hz

% Banda asignada: 2.56kHz (tercera banda - PASABANDA)
fc = 2560; % Hz - Frecuencia central del filtro

%% Cálculo de Frecuencias Límite
% Para encontrar los límites de la banda de paso, calculamos la
% media geométrica de las frecuencias centrales adyacentes según el enunciado.
%
% f_límite = sqrt(f1 * f2)

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
fprintf('Frecuencia de rechazo superior (f2s, -15dB): %.2f Hz\n', f2_15dB);

%% Conversión a Radianes/Segundo (Dominio Analógico)
% Para el diseño analógico, convertimos de Hz a rad/s:
% omega = 2*pi*f

w1 = 2*pi*f1_3dB;   % rad/s - Límite inferior de banda de paso
w2 = 2*pi*f2_3dB;   % rad/s - Límite superior de banda de paso
ws1 = 2*pi*f1_15dB; % rad/s - Frecuencia de rechazo inferior
ws2 = 2*pi*f2_15dB; % rad/s - Frecuencia de rechazo superior

% Frecuencia central del pasabanda (media geométrica en rad/s)
% omega_o = sqrt(omega1 * omega2)
Wo = sqrt(w1*w2);  % rad/s
fc_calculada = Wo/(2*pi);

% Ancho de banda
% BW = omega2 - omega1
BW = w2 - w1;  % rad/s

fprintf('\nFrecuencia central Wo: %.2f rad/s (%.2f Hz)\n', Wo, fc_calculada);
fprintf('Ancho de banda BW: %.2f rad/s\n', BW);

%% Especificaciones de Atenuación
% Según el enunciado del parcial:
Rp = 3;  % Atenuación máxima en la banda de paso (dB)
As = 15; % Atenuación mínima en la banda de rechazo (dB)

fprintf('\nEspecificaciones:\n');
fprintf('  Rp (atenuación en banda de paso): %.2f dB\n', Rp);
fprintf('  As (atenuación en banda de rechazo): %.2f dB\n', As);

%% Cálculo del Orden del Filtro
% Usando buttord() para calcular el orden mínimo requerido.
% Se debe calcular para cada límite y usar el mayor.
%
% La función buttord(Wp, Ws, Rp, As, 's') retorna:
%   - N: orden del filtro
%   - Wc: frecuencia de corte en rad/s

fprintf('\n--- Cálculo del Orden del Filtro ---\n');

% Para el límite inferior (ws1 < w1)
[N_inferior, Wc_inferior] = buttord(w1, ws1, Rp, As, 's');
fprintf('Límite inferior: N = %d, Wc = %.2f rad/s (%.2f Hz)\n', ...
        N_inferior, Wc_inferior, Wc_inferior/(2*pi));

% Para el límite superior (ws2 > w2)
[N_superior, Wc_superior] = buttord(w2, ws2, Rp, As, 's');
fprintf('Límite superior: N = %d, Wc = %.2f rad/s (%.2f Hz)\n', ...
        N_superior, Wc_superior, Wc_superior/(2*pi));

% Seleccionar el orden mayor
N = max(N_inferior, N_superior);
fprintf('\n>>> Orden seleccionado para el filtro: N = %d\n', N);
fprintf('    (El pasabanda tendrá orden 2N = %d)\n', 2*N);

%% Diseño del Filtro Butterworth Pasabajos Prototipo
% Diseñar filtro pasabajos normalizado con orden N y Wc=1
%
% La función butter(N, Wc, 's') diseña un filtro Butterworth analógico:
%   - [b, a] son los coeficientes del numerador y denominador
%   - 's' indica dominio analógico (Laplace)

fprintf('\n--- Diseño del Filtro Prototipo Pasabajos ---\n');
[b_lp, a_lp] = butter(N, 1, 's');
fprintf('Coeficientes del filtro pasabajos prototipo:\n');
fprintf('b = '); disp(b_lp');
fprintf('a = '); disp(a_lp');

%% Transformación a Pasabanda
% Usando lp2bp() para transformar el prototipo pasabajos a pasabanda
%
% La función lp2bp(b, a, Wo, BW) transforma el filtro:
%   - Wo: frecuencia central en rad/s
%   - BW: ancho de banda en rad/s
%   - El orden del pasabanda será 2N

fprintf('\n--- Transformación a Pasabanda ---\n');
[b_bp, a_bp] = lp2bp(b_lp, a_lp, Wo, BW);

fprintf('Coeficientes del filtro pasabanda:\n');
fprintf('Numerador b (orden %d):\n', length(b_bp)-1);
disp(b_bp);
fprintf('Denominador a (orden %d):\n', length(a_bp)-1);
disp(a_bp);

%% Respuesta en Frecuencia del Filtro
% Calcular la respuesta en frecuencia usando freqs()
%
% La función freqs(b, a, w) calcula H(jω) para el filtro analógico

fprintf('\n--- Generando Respuestas en Frecuencia ---\n');

% Vector de frecuencias para análisis (en rad/s)
w_plot = logspace(log10(2*pi*50), log10(2*pi*20000), 1000);

% Calcular respuesta en frecuencia
[H, w_resp] = freqs(b_bp, a_bp, w_plot);

% Convertir a Hz para las gráficas
f_plot = w_resp/(2*pi);

% Magnitud en dB
mag_dB = 20*log10(abs(H));

% Fase en grados
fase_deg = angle(H)*180/pi;

%% Gráfica de Respuesta en Magnitud
figure('Position', [100, 100, 1200, 800]);

% Respuesta en magnitud
subplot(2,1,1);
semilogx(f_plot, mag_dB, 'b', 'LineWidth', 2);
grid on;
hold on;

% Marcar puntos importantes
plot(fc, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Frecuencia central
plot(f1_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2); % -3dB inferior
plot(f2_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2); % -3dB superior
plot(f1_15dB, -15, 'mo', 'MarkerSize', 8, 'LineWidth', 2); % -15dB inferior
plot(f2_15dB, -15, 'mo', 'MarkerSize', 8, 'LineWidth', 2); % -15dB superior

% Líneas de referencia
yline(-3, '--r', 'LineWidth', 1.5);
yline(-15, '--m', 'LineWidth', 1.5);
xline(fc, '--k', 'LineWidth', 1);

xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title(sprintf('Respuesta en Magnitud - Filtro Butterworth Pasabanda (N=%d, fc=%.0f Hz)', N, fc));
legend('Respuesta del filtro', 'Frec. central (0dB)', ...
       'Límites -3dB', '', 'Límites -15dB', '', ...
       '-3dB ref', '-15dB ref', 'fc', 'Location', 'southwest');
xlim([50 20000]);
ylim([-60 5]);

%% Gráfica de Respuesta en Fase
subplot(2,1,2);
semilogx(f_plot, fase_deg, 'r', 'LineWidth', 2);
grid on;
hold on;
xline(fc, '--k', 'LineWidth', 1);

xlabel('Frecuencia (Hz)');
ylabel('Fase (grados)');
title('Respuesta en Fase - Filtro Butterworth Pasabanda');
xlim([50 20000]);

%% Verificación de Especificaciones
% Verificar que el filtro diseñado cumple con las especificaciones
fprintf('\n--- Verificación de Especificaciones ---\n');

% Atenuación en la frecuencia central
idx_fc = find(f_plot >= fc, 1);
atten_fc = mag_dB(idx_fc);
fprintf('Atenuación en fc = %.0f Hz: %.2f dB (esperado: ~0 dB)\n', fc, atten_fc);

% Atenuación en los límites de -3dB
idx_f1 = find(f_plot >= f1_3dB, 1);
atten_f1 = mag_dB(idx_f1);
fprintf('Atenuación en f1 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f1_3dB, atten_f1);

idx_f2 = find(f_plot >= f2_3dB, 1);
atten_f2 = mag_dB(idx_f2);
fprintf('Atenuación en f2 = %.2f Hz: %.2f dB (esperado: -3 dB)\n', f2_3dB, atten_f2);

% Atenuación en las frecuencias de rechazo
idx_fs1 = find(f_plot >= f1_15dB, 1);
atten_fs1 = mag_dB(idx_fs1);
fprintf('Atenuación en fs1 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f1_15dB, atten_fs1);

idx_fs2 = find(f_plot >= f2_15dB, 1);
atten_fs2 = mag_dB(idx_fs2);
fprintf('Atenuación en fs2 = %.0f Hz: %.2f dB (esperado: <= -15 dB)\n', f2_15dB, atten_fs2);

% Verificar cumplimiento
if atten_fs1 <= -15 && atten_fs2 <= -15
    fprintf('\n✓ Las especificaciones de atenuación se cumplen correctamente.\n');
else
    fprintf('\n✗ ADVERTENCIA: Las especificaciones no se cumplen completamente.\n');
end

%% Guardar Resultados
% Guardar coeficientes y parámetros del filtro en archivo .mat
fprintf('\n--- Guardando Resultados ---\n');
save('filtro_butterworth_2560Hz.mat', 'b_bp', 'a_bp', 'N', 'Wo', 'BW', ...
     'fc', 'f1_3dB', 'f2_3dB', 'f1_15dB', 'f2_15dB', 'Rp', 'As');
fprintf('Coeficientes guardados en: filtro_butterworth_2560Hz.mat\n');

% Guardar figura
saveas(gcf, 'respuesta_butterworth_2560Hz.png');
fprintf('Gráfica guardada en: respuesta_butterworth_2560Hz.png\n');

%% Resumen Final
% Resumen de los parámetros del filtro diseñado
fprintf('\n=========================================================================\n');
fprintf('                          RESUMEN DEL DISEÑO\n');
fprintf('=========================================================================\n');
fprintf('Filtro: Butterworth Pasabanda\n');
fprintf('Frecuencia central: %.0f Hz\n', fc);
fprintf('Orden del prototipo pasabajos: N = %d\n', N);
fprintf('Orden del filtro pasabanda: 2N = %d\n', 2*N);
fprintf('Frecuencia central Wo: %.2f rad/s (%.2f Hz)\n', Wo, Wo/(2*pi));
fprintf('Ancho de banda BW: %.2f rad/s (%.2f Hz)\n', BW, BW/(2*pi));
fprintf('Banda de paso (-3dB): %.2f Hz - %.2f Hz\n', f1_3dB, f2_3dB);
fprintf('Frecuencias de rechazo (-15dB): %.0f Hz y %.0f Hz\n', f1_15dB, f2_15dB);
fprintf('=========================================================================\n');

%% Observaciones
% Notas importantes sobre el diseño:
%
% - El filtro pasabanda tiene orden 2N = 4 (prototipo pasabajos N = 2)
% - Las frecuencias de corte Wc calculadas por buttord() son las que
%   garantizan el cumplimiento de las especificaciones
% - Este filtro se usará en el ecualizador de 4 bandas para la frecuencia
%   de 2.56kHz
% - Para el sistema completo, se debe usar el orden mayor de todos los filtros
%
% Funciones MATLAB utilizadas:
% - buttord() - Calcula orden y frecuencia de corte
% - butter() - Diseña filtro Butterworth prototipo
% - lp2bp() - Transforma pasabajos a pasabanda
% - freqs() - Calcula respuesta en frecuencia
