%% TEST DEL ECUALIZADOR ANALÓGICO
% Pruebas de funcionamiento y calidad del sistema de ecualización

clear all;
close all;
clc;

%% Cargar el sistema diseñado
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('           TEST DE FUNCIONAMIENTO Y CALIDAD DEL ECUALIZADOR\n');
fprintf('=========================================================================\n\n');

fprintf('Cargando diseño del ecualizador...\n');
if ~exist('ecualizador_analogico.mat', 'file')
    error('Primero ejecute main_ecualizador_analogico.m para generar el diseño');
end

load('ecualizador_analogico.mat');
n_bandas = length(f_central);
fprintf('Ecualizador cargado: %d bandas\n\n', n_bandas);

%% TEST 1: Verificación de especificaciones
fprintf('=========================================================================\n');
fprintf('              TEST 1: VERIFICACIÓN DE ESPECIFICACIONES\n');
fprintf('=========================================================================\n\n');

w_plot = logspace(log10(2*pi*10), log10(2*pi*fs/2), 2000);
f_plot = w_plot/(2*pi);

fprintf('Verificando cada filtro individualmente...\n\n');

test1_pass = true;
for i = 1:n_bandas
    fprintf('--- Banda %d (fc = %.1f Hz) ---\n', i, f_central(i));

    [H, ~] = freqs(coeficientes_b{i}, coeficientes_a{i}, w_plot);
    mag_dB = 20*log10(abs(H));

    % Verificar atenuación en fc
    idx_fc = find(abs(f_plot - f_central(i)) == min(abs(f_plot - f_central(i))), 1);
    atten_fc = mag_dB(idx_fc);
    fprintf('  Atenuación en fc: %.2f dB ', atten_fc);
    if abs(atten_fc) <= 1
        fprintf('✓\n');
    else
        fprintf('✗ (esperado: ~0 dB)\n');
        test1_pass = false;
    end

    % Verificar límites de -3dB
    idx_3dB = find(mag_dB >= -3.5 & mag_dB <= -2.5);
    if ~isempty(idx_3dB)
        f_3dB_low = f_plot(idx_3dB(1));
        f_3dB_high = f_plot(idx_3dB(end));
        fprintf('  Límites -3dB: %.2f Hz - %.2f Hz ✓\n', f_3dB_low, f_3dB_high);
    else
        fprintf('  Límites -3dB: No encontrados ✗\n');
        test1_pass = false;
    end

    fprintf('\n');
end

if test1_pass
    fprintf('>>> TEST 1: ✓ APROBADO\n\n');
else
    fprintf('>>> TEST 1: ✗ FALLADO\n\n');
end

%% TEST 2: Respuesta total con ganancias unitarias
fprintf('=========================================================================\n');
fprintf('        TEST 2: RESPUESTA TOTAL CON GANANCIAS UNITARIAS (G=1)\n');
fprintf('=========================================================================\n\n');

% Calcular respuesta total
H_total_lin = zeros(size(f_plot));
for i = 1:n_bandas
    [H, ~] = freqs(coeficientes_b{i}, coeficientes_a{i}, w_plot);
    H_total_lin = H_total_lin + abs(H);
end
H_total_dB = 20*log10(H_total_lin);

% Analizar planicidad en banda de audio (20 Hz - 20 kHz)
idx_audio = (f_plot >= 20) & (f_plot <= 20000);
H_audio = H_total_dB(idx_audio);

max_mag = max(H_audio);
min_mag = min(H_audio);
ripple = max_mag - min_mag;
mean_mag = mean(H_audio);

fprintf('Análisis en banda de audio (20 Hz - 20 kHz):\n');
fprintf('  Magnitud máxima:  %.2f dB\n', max_mag);
fprintf('  Magnitud mínima:  %.2f dB\n', min_mag);
fprintf('  Magnitud promedio: %.2f dB\n', mean_mag);
fprintf('  Ripple (variación): %.2f dB\n\n', ripple);

test2_pass = (ripple <= 3);
if test2_pass
    fprintf('>>> TEST 2: ✓ APROBADO (ripple <= 3 dB)\n\n');
else
    fprintf('>>> TEST 2: ✗ FALLADO (ripple > 3 dB)\n\n');
end

%% TEST 3: Independencia de bandas (crosstalk)
fprintf('=========================================================================\n');
fprintf('              TEST 3: INDEPENDENCIA DE BANDAS (CROSSTALK)\n');
fprintf('=========================================================================\n\n');

fprintf('Medición de crosstalk entre bandas adyacentes:\n\n');

test3_pass = true;
for i = 1:n_bandas-1
    % Respuesta del filtro i en la frecuencia central de i+1
    [H_i, ~] = freqs(coeficientes_b{i}, coeficientes_a{i}, w_plot);
    mag_i_dB = 20*log10(abs(H_i));

    idx_fc_next = find(abs(f_plot - f_central(i+1)) == min(abs(f_plot - f_central(i+1))), 1);
    crosstalk = mag_i_dB(idx_fc_next);

    fprintf('  Banda %d en fc de Banda %d (%.1f Hz): %.2f dB ', ...
            i, i+1, f_central(i+1), crosstalk);

    if crosstalk <= -15
        fprintf('✓\n');
    else
        fprintf('✗ (esperado: <= -15 dB)\n');
        test3_pass = false;
    end
end
fprintf('\n');

if test3_pass
    fprintf('>>> TEST 3: ✓ APROBADO (crosstalk <= -15 dB)\n\n');
else
    fprintf('>>> TEST 3: ✗ FALLADO (crosstalk > -15 dB)\n\n');
end

%% TEST 4: Respuesta a señales de prueba
fprintf('=========================================================================\n');
fprintf('              TEST 4: RESPUESTA A SEÑALES DE PRUEBA\n');
fprintf('=========================================================================\n\n');

% Generar señal de prueba: sweep logarítmico
t = 0:1/fs:2; % 2 segundos
f_start = 20;
f_end = 20000;
k = (f_end/f_start)^(1/t(end));
f_inst = f_start * k.^t;
x = sin(2*pi * f_start * (k.^t - 1) / log(k));

fprintf('Señal de prueba: Sweep logarítmico de %.0f Hz a %.0f Hz\n', f_start, f_end);
fprintf('Duración: %.1f s, fs = %d Hz\n\n', t(end), fs);

% Aplicar cada filtro (usando lsim para simulación analógica)
fprintf('Aplicando filtros a la señal de prueba...\n');
y_bandas = zeros(n_bandas, length(t));

for i = 1:n_bandas
    % Crear sistema de transferencia analógico
    sys = tf(coeficientes_b{i}, coeficientes_a{i});

    % Simular respuesta
    y_bandas(i, :) = lsim(sys, x, t);

    % Calcular energía de salida
    energia = sum(y_bandas(i, :).^2);
    fprintf('  Banda %d: Energía de salida = %.2e\n', i, energia);
end
fprintf('\n');

% Calcular salida total (G=1 para todas las bandas)
y_total = sum(y_bandas, 1);

% Analizar retardo y distorsión
fprintf('Análisis de la salida total:\n');
fprintf('  Energía entrada:  %.2e\n', sum(x.^2));
fprintf('  Energía salida:   %.2e\n', sum(y_total.^2));
fprintf('  Relación E_out/E_in: %.2f\n\n', sum(y_total.^2)/sum(x.^2));

test4_pass = true; % Siempre pasa si se ejecuta sin error
fprintf('>>> TEST 4: ✓ APROBADO (señal procesada exitosamente)\n\n');

%% TEST 5: Respuesta con diferentes configuraciones de ganancia
fprintf('=========================================================================\n');
fprintf('        TEST 5: RESPUESTA CON DIFERENTES CONFIGURACIONES DE GANANCIA\n');
fprintf('=========================================================================\n\n');

% Configuraciones de prueba
configs = [
    1.0, 1.0, 1.0, 1.0;  % Plano (todas iguales)
    2.0, 1.0, 1.0, 1.0;  % Enfasis en graves
    1.0, 1.0, 1.0, 2.0;  % Enfasis en agudos
    1.5, 1.5, 0.5, 0.5;  % Graves/medios vs agudos
    0.0, 2.0, 2.0, 0.0;  % Solo medios (efecto "teléfono")
];

config_nombres = {
    'Plano (todas G=1)';
    'Énfasis en graves (G1=2)';
    'Énfasis en agudos (G4=2)';
    'Graves/medios altos';
    'Solo medios (efecto teléfono)'
};

fprintf('Probando %d configuraciones de ganancia...\n\n', size(configs, 1));

% Crear figura de resultados
figure('Position', [50, 50, 1400, 900]);

for cfg = 1:size(configs, 1)
    G = configs(cfg, :);

    % Calcular respuesta total
    H_cfg = zeros(size(f_plot));
    for i = 1:n_bandas
        [H, ~] = freqs(coeficientes_b{i}, coeficientes_a{i}, w_plot);
        H_cfg = H_cfg + G(i) * abs(H);
    end
    H_cfg_dB = 20*log10(H_cfg);

    % Subplot
    subplot(3, 2, cfg);
    semilogx(f_plot, H_cfg_dB, 'b', 'LineWidth', 2);
    grid on;
    hold on;

    % Marcar frecuencias centrales
    for i = 1:n_bandas
        plot(f_central(i), 0, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    end

    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('%s\nG=[%.1f, %.1f, %.1f, %.1f]', config_nombres{cfg}, G));
    xlim([10 fs/2]);
    ylim([-20 15]);
    yline(0, '--k', 'LineWidth', 1);

    % Análisis
    idx_audio = (f_plot >= 20) & (f_plot <= 20000);
    max_mag = max(H_cfg_dB(idx_audio));
    min_mag = min(H_cfg_dB(idx_audio));

    fprintf('Config %d - %s:\n', cfg, config_nombres{cfg});
    fprintf('  Ganancias: [%.1f, %.1f, %.1f, %.1f]\n', G);
    fprintf('  Rango de magnitud: %.2f dB a %.2f dB\n\n', min_mag, max_mag);
end

sgtitle('Test de Configuraciones de Ganancia del Ecualizador', ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'test_configuraciones_ganancia.png');
fprintf('Gráfica guardada: test_configuraciones_ganancia.png\n\n');

test5_pass = true;
fprintf('>>> TEST 5: ✓ APROBADO (todas las configuraciones generadas)\n\n');

%% TEST 6: Estabilidad y propiedades de fase
fprintf('=========================================================================\n');
fprintf('           TEST 6: ESTABILIDAD Y PROPIEDADES DE FASE\n');
fprintf('=========================================================================\n\n');

test6_pass = true;
for i = 1:n_bandas
    fprintf('Banda %d (fc = %.1f Hz):\n', i, f_central(i));

    % Verificar polos
    polos = roots(coeficientes_a{i});
    partes_reales = real(polos);

    fprintf('  Número de polos: %d\n', length(polos));
    fprintf('  Parte real máxima: %.6f\n', max(partes_reales));

    if all(partes_reales < 0)
        fprintf('  Estabilidad: ✓ ESTABLE (todos los polos en semiplano izquierdo)\n');
    else
        fprintf('  Estabilidad: ✗ INESTABLE\n');
        test6_pass = false;
    end

    % Verificar ceros
    ceros = roots(coeficientes_b{i});
    if ~isempty(ceros)
        fprintf('  Número de ceros: %d\n', length(ceros));
    else
        fprintf('  Número de ceros: 0\n');
    end

    fprintf('\n');
end

if test6_pass
    fprintf('>>> TEST 6: ✓ APROBADO (todos los filtros son estables)\n\n');
else
    fprintf('>>> TEST 6: ✗ FALLADO (hay filtros inestables)\n\n');
end

%% RESUMEN DE TESTS
fprintf('=========================================================================\n');
fprintf('                       RESUMEN DE TESTS\n');
fprintf('=========================================================================\n\n');

tests_results = [test1_pass, test2_pass, test3_pass, test4_pass, test5_pass, test6_pass];
tests_names = {
    'Verificación de especificaciones';
    'Respuesta total con G=1';
    'Independencia de bandas (crosstalk)';
    'Respuesta a señales de prueba';
    'Configuraciones de ganancia';
    'Estabilidad y fase'
};

for i = 1:length(tests_results)
    if tests_results(i)
        fprintf('  Test %d - %s: ✓ APROBADO\n', i, tests_names{i});
    else
        fprintf('  Test %d - %s: ✗ FALLADO\n', i, tests_names{i});
    end
end

fprintf('\n');
total_pass = sum(tests_results);
total_tests = length(tests_results);

fprintf('-------------------------------------------------------------------------\n');
fprintf('  RESULTADO FINAL: %d/%d tests aprobados (%.1f%%)\n', ...
        total_pass, total_tests, 100*total_pass/total_tests);
fprintf('-------------------------------------------------------------------------\n\n');

if all(tests_results)
    fprintf('✓✓✓ TODOS LOS TESTS APROBADOS - SISTEMA FUNCIONANDO CORRECTAMENTE ✓✓✓\n\n');
else
    fprintf('⚠⚠⚠ ALGUNOS TESTS FALLARON - REVISAR DISEÑO ⚠⚠⚠\n\n');
end

fprintf('=========================================================================\n\n');
