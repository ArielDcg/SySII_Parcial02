function disenar_filtro_analogico(fc, mostrar_graficas)
%DISENAR_FILTRO_ANALOGICO Diseña filtro analógico Butterworth para ecualizador de 4 bandas
%   fc: Frecuencia central (102.4, 512, 2560, o 12800 Hz)
%   mostrar_graficas: true/false (opcional, default=true)

if nargin < 1
    error('Se requiere el parámetro fc');
end
if nargin < 2
    mostrar_graficas = true;
end

frecuencias_validas = [102.4, 512, 2560, 12800];
if ~ismember(fc, frecuencias_validas)
    error('Frecuencia central inválida. Debe ser: 102.4, 512, 2560, o 12800 Hz');
end

%% Especificaciones del sistema
fprintf('\n=========================================================================\n');
fprintf('      DISEÑO DE FILTRO ANALÓGICO BUTTERWORTH - fc = %.1f Hz\n', fc);
fprintf('=========================================================================\n\n');

f_central = [102.4, 512, 2560, 12800];
fs = 64000;
Rp = 3;   % Atenuación en banda de paso (dB)
As = 15;  % Atenuación en banda de rechazo (dB)

fprintf('Parámetros: fs=%d Hz, Rp=%.1f dB, As=%.1f dB\n\n', fs, Rp, As);

%% Determinar tipo de filtro
fprintf('PASO 1: DETERMINACIÓN DEL TIPO DE FILTRO\n');
fprintf('-------------------------------------------------------------------------\n');

if fc == 102.4
    tipo_filtro = 'PASABAJOS';
    idx_banda = 1;
elseif fc == 12800
    tipo_filtro = 'PASAALTAS';
    idx_banda = 4;
else
    tipo_filtro = 'PASABANDA';
    idx_banda = find(f_central == fc);
end

fprintf('Tipo: %s\n\n', tipo_filtro);

%% Calcular frecuencias límite
fprintf('PASO 2: CÁLCULO DE FRECUENCIAS LÍMITE\n');
fprintf('-------------------------------------------------------------------------\n');

if strcmp(tipo_filtro, 'PASABAJOS')
    % Pasabajos: fc = 102.4 Hz
    fp = sqrt(f_central(1) * f_central(2));           % -3dB
    fs_rechazo = f_central(2);                        % -15dB
    Wp = 2*pi*fp;
    Ws = 2*pi*fs_rechazo;
    f1_3dB = 0;
    f2_3dB = fp;
    f1_15dB = 0;
    f2_15dB = fs_rechazo;

    fprintf('Frecuencia de paso (fp, -3dB): %.2f Hz\n', fp);
    fprintf('Frecuencia de rechazo (fs, -15dB): %.0f Hz\n', fs_rechazo);
    fprintf('Wp = %.2f rad/s, Ws = %.2f rad/s\n\n', Wp, Ws);

elseif strcmp(tipo_filtro, 'PASAALTAS')
    % Pasaaltas: fc = 12800 Hz
    fp = sqrt(f_central(3) * f_central(4));           % -3dB
    fs_rechazo = f_central(3);                        % -15dB
    Wp = 2*pi*fp;
    Ws = 2*pi*fs_rechazo;
    f1_3dB = fp;
    f2_3dB = fs/2;
    f1_15dB = fs_rechazo;
    f2_15dB = fs/2;

    fprintf('Frecuencia de paso (fp, -3dB): %.2f Hz\n', fp);
    fprintf('Frecuencia de rechazo (fs, -15dB): %.0f Hz\n', fs_rechazo);
    fprintf('Wp = %.2f rad/s, Ws = %.2f rad/s\n\n', Wp, Ws);

else
    % Pasabanda: fc = 512 o 2560 Hz
    f1_3dB = sqrt(f_central(idx_banda-1) * f_central(idx_banda));
    f2_3dB = sqrt(f_central(idx_banda) * f_central(idx_banda+1));
    f1_15dB = f_central(idx_banda-1);
    f2_15dB = f_central(idx_banda+1);

    w1 = 2*pi*f1_3dB;
    w2 = 2*pi*f2_3dB;
    ws1 = 2*pi*f1_15dB;
    ws2 = 2*pi*f2_15dB;

    Wo = sqrt(w1*w2);
    BW = w2 - w1;

    fprintf('Banda de paso (-3dB): %.2f Hz - %.2f Hz\n', f1_3dB, f2_3dB);
    fprintf('Frecuencias de rechazo (-15dB): %.1f Hz y %.1f Hz\n', f1_15dB, f2_15dB);
    fprintf('Wo = %.2f rad/s (%.2f Hz), BW = %.2f rad/s\n\n', Wo, Wo/(2*pi), BW);
end

%% Calcular orden del filtro
fprintf('PASO 3: CÁLCULO DEL ORDEN DEL FILTRO\n');
fprintf('-------------------------------------------------------------------------\n');

if strcmp(tipo_filtro, 'PASABANDA')
    [N_inferior, Wc_inferior] = buttord(w1, ws1, Rp, As, 's');
    [N_superior, Wc_superior] = buttord(w2, ws2, Rp, As, 's');
    N = max(N_inferior, N_superior);

    fprintf('Límite inferior: N=%d, Wc=%.2f rad/s\n', N_inferior, Wc_inferior);
    fprintf('Límite superior: N=%d, Wc=%.2f rad/s\n', N_superior, Wc_superior);
    fprintf('Orden seleccionado: N = %d (pasabanda 2N = %d)\n\n', N, 2*N);
else
    [N, Wc] = buttord(Wp, Ws, Rp, As, 's');
    fprintf('Orden: N = %d, Wc = %.2f rad/s (%.2f Hz)\n\n', N, Wc, Wc/(2*pi));
end

%% Diseñar filtro Butterworth
fprintf('PASO 4: DISEÑO DEL FILTRO BUTTERWORTH\n');
fprintf('-------------------------------------------------------------------------\n');

% Prototipo pasabajos normalizado
[b_lp, a_lp] = butter(N, 1, 's');
fprintf('Filtro prototipo pasabajos (N=%d, Wc=1 rad/s)\n', N);

% Transformación al tipo deseado
if strcmp(tipo_filtro, 'PASABANDA')
    [b, a] = lp2bp(b_lp, a_lp, Wo, BW);
    fprintf('Transformación a pasabanda: lp2bp(b_lp, a_lp, %.2f, %.2f)\n', Wo, BW);
elseif strcmp(tipo_filtro, 'PASAALTAS')
    [b, a] = lp2hp(b_lp, a_lp, Wc);
    fprintf('Transformación a pasaaltas: lp2hp(b_lp, a_lp, %.2f)\n', Wc);
else
    [b, a] = butter(N, Wc, 's');
    fprintf('Desnormalización: butter(%d, %.2f, ''s'')\n', N, Wc);
end

fprintf('Orden del filtro final: %d\n\n', length(a)-1);

%% Calcular respuesta en frecuencia
fprintf('PASO 5: RESPUESTA EN FRECUENCIA\n');
fprintf('-------------------------------------------------------------------------\n');

w_plot = logspace(log10(2*pi*10), log10(2*pi*fs/2), 2000);
[H, w_resp] = freqs(b, a, w_plot);
f_plot = w_resp/(2*pi);
mag_dB = 20*log10(abs(H));
fase_deg = angle(H)*180/pi;

fprintf('Respuesta calculada: 2000 puntos, 10 Hz a %.0f Hz\n\n', fs/2);

%% Verificar especificaciones
fprintf('PASO 6: VERIFICACIÓN DE ESPECIFICACIONES\n');
fprintf('-------------------------------------------------------------------------\n');

find_idx = @(f_vec, f_target) find(abs(f_vec - f_target) == min(abs(f_vec - f_target)), 1);

% Atenuación en fc
idx_fc = find_idx(f_plot, fc);
atten_fc = mag_dB(idx_fc);
fprintf('Atenuación en fc=%.1f Hz: %.2f dB ', fc, atten_fc);
if abs(atten_fc) <= 1
    fprintf('✓\n');
else
    fprintf('✗\n');
end

cumple = abs(atten_fc) <= 1;

% Verificar límites según tipo
if strcmp(tipo_filtro, 'PASABANDA')
    idx_f1 = find_idx(f_plot, f1_3dB);
    atten_f1 = mag_dB(idx_f1);
    if abs(atten_f1 + 3) <= 0.5
        fprintf('Atenuación en f1=%.2f Hz: %.2f dB ✓\n', f1_3dB, atten_f1);
    else
        fprintf('Atenuación en f1=%.2f Hz: %.2f dB ✗\n', f1_3dB, atten_f1);
    end

    idx_f2 = find_idx(f_plot, f2_3dB);
    atten_f2 = mag_dB(idx_f2);
    if abs(atten_f2 + 3) <= 0.5
        fprintf('Atenuación en f2=%.2f Hz: %.2f dB ✓\n', f2_3dB, atten_f2);
    else
        fprintf('Atenuación en f2=%.2f Hz: %.2f dB ✗\n', f2_3dB, atten_f2);
    end

    idx_fs1 = find_idx(f_plot, f1_15dB);
    atten_fs1 = mag_dB(idx_fs1);
    if atten_fs1 <= -15
        fprintf('Atenuación en fs1=%.1f Hz: %.2f dB ✓\n', f1_15dB, atten_fs1);
    else
        fprintf('Atenuación en fs1=%.1f Hz: %.2f dB ✗\n', f1_15dB, atten_fs1);
    end

    idx_fs2 = find_idx(f_plot, f2_15dB);
    atten_fs2 = mag_dB(idx_fs2);
    if atten_fs2 <= -15
        fprintf('Atenuación en fs2=%.1f Hz: %.2f dB ✓\n', f2_15dB, atten_fs2);
    else
        fprintf('Atenuación en fs2=%.1f Hz: %.2f dB ✗\n', f2_15dB, atten_fs2);
    end

    cumple = cumple && (abs(atten_f1 + 3) <= 0.5) && (abs(atten_f2 + 3) <= 0.5) && ...
             (atten_fs1 <= -15) && (atten_fs2 <= -15);

elseif strcmp(tipo_filtro, 'PASABAJOS')
    idx_fp = find_idx(f_plot, f2_3dB);
    atten_fp = mag_dB(idx_fp);
    if abs(atten_fp + 3) <= 0.5
        fprintf('Atenuación en fp=%.2f Hz: %.2f dB ✓\n', f2_3dB, atten_fp);
    else
        fprintf('Atenuación en fp=%.2f Hz: %.2f dB ✗\n', f2_3dB, atten_fp);
    end

    idx_fs = find_idx(f_plot, f2_15dB);
    atten_fs = mag_dB(idx_fs);
    if atten_fs <= -15
        fprintf('Atenuación en fs=%.0f Hz: %.2f dB ✓\n', f2_15dB, atten_fs);
    else
        fprintf('Atenuación en fs=%.0f Hz: %.2f dB ✗\n', f2_15dB, atten_fs);
    end

    cumple = cumple && (abs(atten_fp + 3) <= 0.5) && (atten_fs <= -15);

else % PASAALTAS
    idx_fp = find_idx(f_plot, f1_3dB);
    atten_fp = mag_dB(idx_fp);
    if abs(atten_fp + 3) <= 0.5
        fprintf('Atenuación en fp=%.2f Hz: %.2f dB ✓\n', f1_3dB, atten_fp);
    else
        fprintf('Atenuación en fp=%.2f Hz: %.2f dB ✗\n', f1_3dB, atten_fp);
    end

    idx_fs = find_idx(f_plot, f1_15dB);
    atten_fs = mag_dB(idx_fs);
    if atten_fs <= -15
        fprintf('Atenuación en fs=%.0f Hz: %.2f dB ✓\n', f1_15dB, atten_fs);
    else
        fprintf('Atenuación en fs=%.0f Hz: %.2f dB ✗\n', f1_15dB, atten_fs);
    end

    cumple = cumple && (abs(atten_fp + 3) <= 0.5) && (atten_fs <= -15);
end

fprintf('\n');
if cumple
    fprintf('✓ Todas las especificaciones se cumplen\n\n');
else
    fprintf('✗ Algunas especificaciones no se cumplen\n\n');
end

%% Generar gráficas
if mostrar_graficas
    fprintf('PASO 7: GENERACIÓN DE GRÁFICAS\n');
    fprintf('-------------------------------------------------------------------------\n');

    fig = figure('Position', [100, 100, 1200, 800]);
    fig.Name = sprintf('Filtro %s - fc = %.1f Hz', tipo_filtro, fc);

    % Magnitud
    subplot(2,1,1);
    semilogx(f_plot, mag_dB, 'b', 'LineWidth', 2);
    grid on;
    hold on;
    plot(fc, atten_fc, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    if strcmp(tipo_filtro, 'PASABANDA')
        plot(f1_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2);
        plot(f2_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2);
        plot(f1_15dB, atten_fs1, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
        plot(f2_15dB, atten_fs2, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
    elseif strcmp(tipo_filtro, 'PASABAJOS')
        plot(f2_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2);
        plot(f2_15dB, atten_fs, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
    else
        plot(f1_3dB, -3, 'go', 'MarkerSize', 8, 'LineWidth', 2);
        plot(f1_15dB, atten_fs, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
    end

    yline(-3, '--r', 'LineWidth', 1.5);
    yline(-15, '--m', 'LineWidth', 1.5);
    xline(fc, '--k', 'LineWidth', 1);

    xlabel('Frecuencia (Hz)');
    ylabel('Magnitud (dB)');
    title(sprintf('Magnitud - Butterworth %s (N=%d, fc=%.1f Hz)', tipo_filtro, length(a)-1, fc));
    xlim([10 fs/2]);
    ylim([-80 5]);
    legend('Location', 'southwest');

    % Fase
    subplot(2,1,2);
    semilogx(f_plot, fase_deg, 'r', 'LineWidth', 2);
    grid on;
    hold on;
    xline(fc, '--k', 'LineWidth', 1);
    xlabel('Frecuencia (Hz)');
    ylabel('Fase (grados)');
    title(sprintf('Fase - Butterworth %s', tipo_filtro));
    xlim([10 fs/2]);

    fprintf('Gráficas generadas\n\n');
end

%% Guardar resultados
fprintf('PASO 8: GUARDAR RESULTADOS\n');
fprintf('-------------------------------------------------------------------------\n');

if fc == 102.4
    nombre_base = 'filtro_butterworth_102p4Hz';
else
    nombre_base = sprintf('filtro_butterworth_%.0fHz', fc);
end

nombre_mat = [nombre_base '.mat'];
if strcmp(tipo_filtro, 'PASABANDA')
    save(nombre_mat, 'b', 'a', 'N', 'Wo', 'BW', 'fc', 'tipo_filtro', ...
         'f1_3dB', 'f2_3dB', 'f1_15dB', 'f2_15dB', 'Rp', 'As');
else
    save(nombre_mat, 'b', 'a', 'N', 'fc', 'tipo_filtro', 'Rp', 'As');
end
fprintf('Coeficientes: %s\n', nombre_mat);

if mostrar_graficas
    nombre_png = [nombre_base '_respuesta.png'];
    saveas(gcf, nombre_png);
    fprintf('Gráfica: %s\n', nombre_png);
end

fprintf('\n=========================================================================\n');
if cumple
    fprintf('RESUMEN: %s | fc=%.1f Hz | Orden=%d | Estado: ✓ OK\n', ...
            tipo_filtro, fc, length(a)-1);
else
    fprintf('RESUMEN: %s | fc=%.1f Hz | Orden=%d | Estado: ✗ ERROR\n', ...
            tipo_filtro, fc, length(a)-1);
end
fprintf('=========================================================================\n\n');

end
