%% Ecualizador de Audio de 4 bandas - Sistema FIR
% Segundo Parcial - Señales y Sistemas II
% Sistema completo usando filtros FIR

clear all;
close all;
clc;

%% Cargar filtros diseñados
fprintf('=== ECUALIZADOR DE AUDIO - SISTEMA FIR ===\n\n');
fprintf('Cargando filtros FIR diseñados...\n');

if ~exist('filtros_FIR_diseñados.mat', 'file')
    error('No se encontró el archivo filtros_FIR_diseñados.mat. Ejecute primero diseno_filtros_FIR.m');
end

load('filtros_FIR_diseñados.mat');
fprintf('Filtros cargados exitosamente.\n\n');

%% Configuración de ganancias
% Las ganancias pueden tomar valores entre 0 y 2
% Vector de ganancias [G1, G2, G3, G4] para las 4 bandas

% Ejemplo: [0.25, 0.5, 1, 2] (especificado en el enunciado)
ganancias = [0.25, 0.5, 1, 2];

fprintf('--- Configuración de ganancias ---\n');
fprintf('G1 (%.1f Hz):   %.2f\n', f_central(1), ganancias(1));
fprintf('G2 (%.0f Hz):  %.2f\n', f_central(2), ganancias(2));
fprintf('G3 (%.2f kHz): %.2f\n', f_central(3)/1000, ganancias(3));
fprintf('G4 (%.1f kHz): %.2f\n\n', f_central(4)/1000, ganancias(4));

%% Generar o cargar señal de entrada
% Opción 1: Generar señal de prueba con múltiples frecuencias
% Opción 2: Cargar archivo de audio

opcion_entrada = 1; % 1: Señal de prueba, 2: Archivo de audio

if opcion_entrada == 1
    fprintf('Generando señal de prueba...\n');

    % Duración de la señal
    duracion = 2; % segundos
    t = 0:1/fs:(duracion - 1/fs);
    N_samples = length(t);

    % Crear señal de prueba con componentes en cada banda
    % Banda 1: 102.4 Hz
    f1 = 102.4;
    A1 = 1.0;

    % Banda 2: 512 Hz
    f2 = 512;
    A2 = 0.8;

    % Banda 3: 2560 Hz
    f3 = 2560;
    A3 = 0.6;

    % Banda 4: 12800 Hz
    f4 = 12800;
    A4 = 0.4;

    % Señal compuesta
    x = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t) + ...
        A3*sin(2*pi*f3*t) + A4*sin(2*pi*f4*t);

    % Normalizar
    x = x / max(abs(x));

    fprintf('Señal de prueba generada: %.1f segundos, %d muestras\n', duracion, N_samples);
    fprintf('Componentes de frecuencia: %.1f Hz, %.0f Hz, %.2f kHz, %.1f kHz\n\n', ...
        f1, f2, f3/1000, f4/1000);
else
    fprintf('Buscando archivo de audio...\n');
    % Intentar cargar un archivo de audio
    % Nota: Ajustar la ruta según sea necesario
    [x, fs_audio] = audioread('audio_entrada.wav');

    % Si la frecuencia de muestreo es diferente, remuestrear
    if fs_audio ~= fs
        fprintf('Remuestreando de %d Hz a %d Hz...\n', fs_audio, fs);
        x = resample(x, fs, fs_audio);
    end

    % Si es estéreo, convertir a mono
    if size(x, 2) > 1
        x = mean(x, 2);
    end

    fprintf('Archivo de audio cargado: %d muestras\n\n', length(x));
end

%% Aplicar filtros FIR a la señal de entrada
fprintf('--- Aplicando filtros FIR ---\n');

% Filtro 1: Pasabajos (102.4 Hz)
fprintf('Filtrando banda 1 (%.1f Hz)...\n', f_central(1));
y1 = filter(b1, 1, x);
y1 = ganancias(1) * y1;

% Filtro 2: Pasabanda (512 Hz)
fprintf('Filtrando banda 2 (%.0f Hz)...\n', f_central(2));
y2 = filter(b2, 1, x);
y2 = ganancias(2) * y2;

% Filtro 3: Pasabanda (2560 Hz)
fprintf('Filtrando banda 3 (%.2f kHz)...\n', f_central(3)/1000);
y3 = filter(b3, 1, x);
y3 = ganancias(3) * y3;

% Filtro 4: Pasaaltos (12800 Hz)
fprintf('Filtrando banda 4 (%.1f kHz)...\n', f_central(4)/1000);
y4 = filter(b4, 1, x);
y4 = ganancias(4) * y4;

fprintf('Filtrado completado.\n\n');

%% Combinar señales filtradas (suma)
fprintf('Combinando señales filtradas...\n');
y_total = y1 + y2 + y3 + y4;

% Normalizar salida para evitar saturación
max_val = max(abs(y_total));
if max_val > 1
    fprintf('Normalizando salida (factor: %.2f)...\n', max_val);
    y_total = y_total / max_val;
end

fprintf('Ecualizador aplicado exitosamente.\n\n');

%% Análisis espectral
fprintf('Calculando espectros...\n');

% Parámetros para FFT
Nfft = 2^nextpow2(length(x));
f_fft = (0:Nfft-1) * (fs/Nfft);
f_fft = f_fft(1:Nfft/2+1); % Solo frecuencias positivas

% FFT de la señal de entrada
X = fft(x, Nfft);
X_mag = abs(X(1:Nfft/2+1));
X_mag_dB = 20*log10(X_mag + eps);

% FFT de las salidas de cada filtro
Y1 = fft(y1, Nfft);
Y1_mag = abs(Y1(1:Nfft/2+1));
Y1_mag_dB = 20*log10(Y1_mag + eps);

Y2 = fft(y2, Nfft);
Y2_mag = abs(Y2(1:Nfft/2+1));
Y2_mag_dB = 20*log10(Y2_mag + eps);

Y3 = fft(y3, Nfft);
Y3_mag = abs(Y3(1:Nfft/2+1));
Y3_mag_dB = 20*log10(Y3_mag + eps);

Y4 = fft(y4, Nfft);
Y4_mag = abs(Y4(1:Nfft/2+1));
Y4_mag_dB = 20*log10(Y4_mag + eps);

% FFT de la salida total
Y_total = fft(y_total, Nfft);
Y_total_mag = abs(Y_total(1:Nfft/2+1));
Y_total_mag_dB = 20*log10(Y_total_mag + eps);

fprintf('Espectros calculados.\n\n');

%% Gráficas - Espectros de entrada y salida de cada filtro
figure('Name', 'Espectros - Entrada y Salida de cada Filtro FIR', ...
       'Position', [100 100 1200 900]);

% Filtro 1
subplot(4,2,1);
plot(t(1:1000), x(1:1000), 'b', 'LineWidth', 1);
grid on;
title(sprintf('Filtro 1 (%.1f Hz) - Entrada', f_central(1)));
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(4,2,2);
semilogx(f_fft, X_mag_dB, 'b', 'LineWidth', 1);
hold on;
semilogx(f_fft, Y1_mag_dB, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 1 (%.1f Hz) - Espectro', f_central(1)));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
ylim([min(X_mag_dB)-10 max(X_mag_dB)+10]);
legend('Entrada', 'Salida');

% Filtro 2
subplot(4,2,3);
plot(t(1:1000), x(1:1000), 'b', 'LineWidth', 1);
grid on;
title(sprintf('Filtro 2 (%.0f Hz) - Entrada', f_central(2)));
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(4,2,4);
semilogx(f_fft, X_mag_dB, 'b', 'LineWidth', 1);
hold on;
semilogx(f_fft, Y2_mag_dB, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 2 (%.0f Hz) - Espectro', f_central(2)));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
ylim([min(X_mag_dB)-10 max(X_mag_dB)+10]);
legend('Entrada', 'Salida');

% Filtro 3
subplot(4,2,5);
plot(t(1:1000), x(1:1000), 'b', 'LineWidth', 1);
grid on;
title(sprintf('Filtro 3 (%.2f kHz) - Entrada', f_central(3)/1000));
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(4,2,6);
semilogx(f_fft, X_mag_dB, 'b', 'LineWidth', 1);
hold on;
semilogx(f_fft, Y3_mag_dB, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 3 (%.2f kHz) - Espectro', f_central(3)/1000));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
ylim([min(X_mag_dB)-10 max(X_mag_dB)+10]);
legend('Entrada', 'Salida');

% Filtro 4
subplot(4,2,7);
plot(t(1:1000), x(1:1000), 'b', 'LineWidth', 1);
grid on;
title(sprintf('Filtro 4 (%.1f kHz) - Entrada', f_central(4)/1000));
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(4,2,8);
semilogx(f_fft, X_mag_dB, 'b', 'LineWidth', 1);
hold on;
semilogx(f_fft, Y4_mag_dB, 'r', 'LineWidth', 1.5);
grid on;
title(sprintf('Filtro 4 (%.1f kHz) - Espectro', f_central(4)/1000));
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
ylim([min(X_mag_dB)-10 max(X_mag_dB)+10]);
legend('Entrada', 'Salida');

saveas(gcf, 'FIR_Ecualizador_Espectros_Filtros.png');

%% Gráficas - Entrada y salida del sistema completo
figure('Name', 'Sistema Completo - Ecualizador FIR', 'Position', [200 200 1200 700]);

% Señales en el tiempo
subplot(3,1,1);
plot(t, x, 'b', 'LineWidth', 1);
grid on;
title('Señal de Entrada');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 min(0.1, max(t))]);

subplot(3,1,2);
plot(t, y_total, 'r', 'LineWidth', 1);
grid on;
title(sprintf('Señal de Salida Ecualizada - Ganancias: [%.2f, %.2f, %.2f, %.2f]', ...
      ganancias(1), ganancias(2), ganancias(3), ganancias(4)));
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 min(0.1, max(t))]);

% Espectros
subplot(3,1,3);
semilogx(f_fft, X_mag_dB, 'b', 'LineWidth', 1.5);
hold on;
semilogx(f_fft, Y_total_mag_dB, 'r', 'LineWidth', 1.5);
grid on;
xline(f_central(1), 'k--', sprintf('%.1f Hz', f_central(1)));
xline(f_central(2), 'k--', sprintf('%.0f Hz', f_central(2)));
xline(f_central(3), 'k--', sprintf('%.2f kHz', f_central(3)/1000));
xline(f_central(4), 'k--', sprintf('%.1f kHz', f_central(4)/1000));
title('Espectros de Magnitud - Entrada vs Salida');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
legend('Entrada', 'Salida Ecualizada', 'Location', 'southwest');

saveas(gcf, 'FIR_Ecualizador_Sistema_Completo.png');

%% Gráfica - Todas las bandas filtradas
figure('Name', 'Bandas Filtradas Individuales - FIR', 'Position', [250 250 1200 600]);

subplot(2,1,1);
plot(t, y1, 'b', 'LineWidth', 1); hold on;
plot(t, y2, 'r', 'LineWidth', 1);
plot(t, y3, 'g', 'LineWidth', 1);
plot(t, y4, 'm', 'LineWidth', 1);
grid on;
title('Salidas de los Filtros (en el tiempo)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
xlim([0 min(0.1, max(t))]);
legend(sprintf('Banda 1 (G=%.2f)', ganancias(1)), ...
       sprintf('Banda 2 (G=%.2f)', ganancias(2)), ...
       sprintf('Banda 3 (G=%.2f)', ganancias(3)), ...
       sprintf('Banda 4 (G=%.2f)', ganancias(4)));

subplot(2,1,2);
semilogx(f_fft, Y1_mag_dB, 'b', 'LineWidth', 1.5); hold on;
semilogx(f_fft, Y2_mag_dB, 'r', 'LineWidth', 1.5);
semilogx(f_fft, Y3_mag_dB, 'g', 'LineWidth', 1.5);
semilogx(f_fft, Y4_mag_dB, 'm', 'LineWidth', 1.5);
grid on;
title('Espectros de las Bandas Filtradas');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
xlim([10 fs/2]);
legend('Banda 1 (102.4 Hz)', 'Banda 2 (512 Hz)', ...
       'Banda 3 (2.56 kHz)', 'Banda 4 (12.8 kHz)');

saveas(gcf, 'FIR_Ecualizador_Bandas_Individuales.png');

%% Guardar resultados
fprintf('Guardando resultados...\n');
save('ecualizador_FIR_resultados.mat', 'x', 'y1', 'y2', 'y3', 'y4', ...
     'y_total', 'ganancias', 'fs', 'f_central');

% Guardar señal de salida como audio (opcional)
if opcion_entrada == 1
    % Normalizar a rango [-1, 1] para audio
    y_audio = y_total / max(abs(y_total));
    audiowrite('salida_ecualizada_FIR.wav', y_audio, fs);
    fprintf('Audio de salida guardado: salida_ecualizada_FIR.wav\n');
end

fprintf('\n=== ECUALIZADOR FIR COMPLETADO ===\n');
fprintf('Archivos generados:\n');
fprintf('  - ecualizador_FIR_resultados.mat\n');
fprintf('  - FIR_Ecualizador_Espectros_Filtros.png\n');
fprintf('  - FIR_Ecualizador_Sistema_Completo.png\n');
fprintf('  - FIR_Ecualizador_Bandas_Individuales.png\n');
if opcion_entrada == 1
    fprintf('  - salida_ecualizada_FIR.wav\n');
end
fprintf('\n');

%% Estadísticas del sistema
fprintf('--- Estadísticas del sistema ---\n');
fprintf('Número total de coeficientes: %d\n', length(b1) + length(b2) + length(b3) + length(b4));
fprintf('Retardo total del sistema (máximo): %d muestras (%.3f ms)\n', ...
        max([N1, N2, N3, N4])/2, max([N1, N2, N3, N4])/(2*fs)*1000);
fprintf('Fase lineal: Sí (filtros FIR con simetría)\n');
fprintf('Estabilidad: Garantizada (filtros FIR)\n\n');
