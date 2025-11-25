# Filtro Butterworth Pasabanda para 2.56kHz

## Descripción
Este script diseña un filtro Butterworth pasabanda analógico para la frecuencia de 2.56kHz como parte del ecualizador de audio de 4 bandas.

## Archivos
- `butterworth_pasabanda_2560Hz.m`: Script principal de diseño
- `filtro_butterworth_2560Hz.mat`: Coeficientes del filtro (generado al ejecutar)
- `respuesta_butterworth_2560Hz.png`: Gráfica de respuesta en frecuencia (generada al ejecutar)

## Cómo usar

### Ejecutar en MATLAB
```matlab
% Simplemente ejecuta el script
butterworth_pasabanda_2560Hz
```

### Lo que hace el script:
1. **Define las especificaciones** según el enunciado del parcial
2. **Calcula las frecuencias límite** usando la media geométrica
3. **Calcula el orden del filtro** usando `buttord()`
4. **Diseña el filtro prototipo pasabajos** con `butter()`
5. **Transforma a pasabanda** usando `lp2bp()`
6. **Genera gráficas** de magnitud y fase
7. **Verifica** que se cumplan las especificaciones
8. **Guarda los resultados** en archivos .mat y .png

## Especificaciones del Filtro

### Frecuencias
- **Frecuencia central (fc)**: 2560 Hz
- **Banda de paso (-3dB)**: ~1145 Hz - ~5724 Hz
- **Banda de rechazo (-15dB)**: 512 Hz y 12800 Hz

### Requisitos
- Atenuación en banda de paso: **3 dB** (en los límites)
- Atenuación en banda de rechazo: **15 dB** (en frecuencias centrales adyacentes)
- Ganancia a frecuencia central: **0 dB**

## Resultados Esperados

### Orden del Filtro
Basado en tus cálculos:
- **Orden del prototipo pasabajos**: N = 2
- **Orden del filtro pasabanda**: 2N = 4

### Frecuencias de Corte
Según tus valores:
- **Wc1 = 3168.09 Hz** (límite inferior)
- **Wc2 = 4353.00 Hz** (límite superior)

Estos valores corresponden a las frecuencias de corte del filtro pasabanda donde la atenuación es -3dB.

## Fórmulas Utilizadas

### Media Geométrica (límites de banda de paso)
```
f1 = √(512 × 2560) ≈ 1145 Hz
f2 = √(2560 × 12800) ≈ 5724 Hz
```

### Frecuencia Central (Wo)
```
Wo = √(w1 × w2)  [rad/s]
donde w1 = 2πf1 y w2 = 2πf2
```

### Ancho de Banda (BW)
```
BW = w2 - w1  [rad/s]
```

### Orden del Filtro
```
N ≥ (1/2) × log₁₀[(10^(Rp/10) - 1)/(10^(As/10) - 1)] / log₁₀(w1/w2)
```

## Funciones MATLAB Utilizadas

1. **`buttord(W1, W2, Rp, As, 's')`**
   - Calcula el orden N y frecuencia de corte Wc del filtro

2. **`butter(N, Wc, 's')`**
   - Diseña un filtro Butterworth analógico de orden N

3. **`lp2bp(b, a, Wo, BW)`**
   - Transforma un filtro pasabajos a pasabanda
   - Wo: frecuencia central [rad/s]
   - BW: ancho de banda [rad/s]

4. **`freqs(b, a, w)`**
   - Calcula la respuesta en frecuencia del filtro analógico

## Notas Importantes

1. **Dominio Analógico**: Este filtro está diseñado en el dominio analógico (parámetro 's' en las funciones)

2. **Orden del Pasabanda**: El filtro pasabanda tiene orden 2N, donde N es el orden del prototipo pasabajos

3. **Verificación**: El script verifica automáticamente que se cumplan las especificaciones de atenuación

4. **Extensión**: Este diseño se puede extender para los otros filtros del ecualizador (102.4Hz, 512Hz, 12.8kHz)

## Para el Informe

El script genera toda la información necesaria para el informe:
- ✓ Procedimientos de cálculo
- ✓ Orden del filtro
- ✓ Frecuencias de corte
- ✓ Coeficientes del numerador y denominador
- ✓ Gráficas de magnitud y fase
- ✓ Verificación de especificaciones

## Siguiente Paso

Para completar el punto 1 del parcial (filtros analógicos de Butterworth), necesitas:
1. Diseñar los otros 3 filtros (102.4Hz pasabajos, 512Hz pasabanda, 12.8kHz pasaaltos)
2. Usar el **orden mayor** de todos los filtros para diseñarlos todos
3. Implementar el sistema completo del ecualizador

---
**Autor**: Diseño para frecuencia 2.56kHz
**Materia**: Señales y Sistemas II - Segundo Parcial
**Fecha**: 2025
