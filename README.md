# 🌫️ SEVA-WRF — Sistema de Evaluación y Visualización Avanzada del Pronóstico WRF-Chem

[![bash](https://img.shields.io/badge/bash-%E2%89%A54.0-blue?logo=gnu-bash)](#requisitos-del-sistema)
[![python](https://img.shields.io/badge/python-%E2%89%A53.8-blue?logo=python)](#dependencias)
[![license](https://img.shields.io/badge/license-MIT-green)](#licencia)
[![version](https://img.shields.io/badge/versi%C3%B3n-1.0.0-orange)](#changelog)

Suite de herramientas Python para el **análisis estadístico avanzado**, **visualización del desempeño** y **pronóstico ponderado** del sistema de calidad del aire WRF-Chem, comparado contra observaciones horarias de la red **SINAICA/INECC** en ocho zonas metropolitanas del centro de México.

SEVA-WRF complementa el pipeline operativo [ddsinaica](https://github.com/JoseAgustin/ddsinaica) con módulos de diagnóstico científico profundo: análisis por percentiles locales, curvas de habilidad multi-umbral, diagramas de Taylor y de Roebber, mapas de calor temporales y combinación ponderada de horizontes de pronóstico.

---

## Tabla de contenidos

- [Descripción](#descripción)
- [Arquitectura](#arquitectura)
- [Requisitos e instalación](#requisitos-e-instalación)
- [Scripts Python](#scripts-python)
  - [taylor_mensual.py](#1-taylor_mensualpy)
  - [heatmap_desempeno.py](#2-heatmap_desempenopy)
  - [roebber_desempeno.py](#3-roebber_desempenopy)
  - [informe_dicotomico.py](#4-informe_dicotomicopy)
  - [informe_graficas.py](#5-informe_graficaspy)
  - [percentiles_obs.py](#6-percentiles_obspy)
  - [curva_habilidad.py](#7-curva_habilidadpy)
  - [analisis_espacial_percentil.py](#8-analisis_espacial_percentilpy)
  - [pronostico_ponderado.py](#9-pronostico_ponderadopy)
- [Ciudades y contaminantes](#ciudades-y-contaminantes)
- [Umbrales normativos](#umbrales-normativos)
- [Temporadas climáticas](#temporadas-climáticas)
- [Flujo de trabajo recomendado](#flujo-de-trabajo-recomendado)
- [Changelog](#changelog)
- [Licencia](#licencia)

---

## Descripción

SEVA-WRF implementa nueve módulos de análisis organizados en tres niveles:

| Nivel | Módulo | Función principal |
|-------|--------|-------------------|
| **Diagnóstico continuo** | `taylor_mensual.py` | Diagramas de Taylor mensuales (σ, R, CRMSE normalizados) |
| **Diagnóstico dicotómico** | `informe_dicotomico.py` | Informe Word con POD, FAR, CSI, TSS, PC, BIAS por ciudad × mes |
| **Visualización temporal** | `heatmap_desempeno.py` | Mapas de calor del desempeño por ciudad y mes |
| **Visualización dicotómica** | `roebber_desempeno.py` | Diagrama de Roebber (trayectoria POD vs SR) |
| **Informe integrado** | `informe_graficas.py` | Documento Word con todas las visualizaciones y sus descripciones |
| **Análisis de percentiles** | `percentiles_obs.py` | Percentiles empíricos locales y comparación modelo vs obs |
| **Curva de habilidad** | `curva_habilidad.py` | Sensibilidad del pronóstico a umbrales alternativos |
| **Análisis espacial** | `analisis_espacial_percentil.py` | Evaluación dicotómica con umbrales adaptativos locales |
| **Pronóstico ponderado** | `pronostico_ponderado.py` | Combinación de horizontes +24h/+48h/+72h por RMSE inverso |

---

## Arquitectura

```
combinado/ajustados/
└── eval_<CONT>_<Ciudad>_YYYY-MM-DD.csv
    Fecha, Ciudad, max_obs, mod_dia1, mod_dia2, mod_dia3
         │
         ├─► taylor_mensual.py          → resultados_taylor/
         │     taylor_YYYY_MM.png
         │     estadisticas_taylor.csv
         │
         ├─► heatmap_desempeno.py       → resultados_heatmap/
         │     heatmap_<CONT>_<MET>_<HOR>_<RES>_<CAT>.png
         │
         ├─► roebber_desempeno.py       → resultados_roebber/
         │     roebber_<CONT>_<AGRUP>_<HOR>_<CAT>.png
         │
         ├─► informe_dicotomico.py      → informes_dicotomicos/
         │     informe_dicotomico_<CAT>.docx
         │     dicotomico_stats.csv
         │
         ├─► percentiles_obs.py         → resultados_percentiles/
         │     percentiles_<CONT>_P<N>_heatmap_obs_<HOR>.png
         │     percentiles_<CONT>_P<N>_heatmap_sesgo_<HOR>.png
         │     percentiles_<CONT>_QQ_<HOR>.png
         │     percentiles_<CONT>_serie_<Ciudad>_<HOR>.png
         │     tabla_percentiles_<CONT>.csv
         │
         ├─► curva_habilidad.py         → resultados_curva_habilidad/
         │     curva_habilidad_<CONT>_<HOR>_<AGRUP>.png
         │     tabla_barrido_umbrales.csv
         │
         ├─► analisis_espacial_percentil.py → resultados_espacial_percentil/
         │     espacial_<CONT>_<HOR>_<MODO>_P<N>.png
         │     espacial_<CONT>_<HOR>_umbral_P<N>_<MODO>.png
         │     espacial_<CONT>_<HOR>_delta_P<N>_<MODO>.png
         │     espacial_<CONT>_<HOR>_comparacion_P<N>_<MODO>.png
         │     tabla_espacial_percentil_<CONT>.csv
         │
         ├─► pronostico_ponderado.py    → resultados_ponderado/
         │     pronostico_ponderado_HOY.csv
         │     pronostico_ponderado_MANANA.csv
         │     pesos_historicos.csv
         │     serie_ponderada_historica.csv
         │     pesos_por_temporada.png
         │     pronostico_ponderado_<CONT>_serie.png
         │     pronostico_ponderado_scatter.png
         │
         └─► informe_graficas.py        → informe_visualizaciones.docx
               (recolecta todos los PNG y genera documento Word integrado)
```

---

## Requisitos e instalación

### Python ≥ 3.8

```bash
pip install pandas numpy matplotlib scipy python-docx Pillow
```

| Paquete | Uso |
|---------|-----|
| `pandas` | Lectura y manipulación de CSV |
| `numpy` | Cálculo numérico y estadístico |
| `matplotlib` | Generación de todas las gráficas |
| `scipy` | Correlación de Pearson, regresión, filtrado |
| `python-docx` | Generación de documentos Word (.docx) |
| `Pillow` | Inserción de imágenes en documentos Word |

### Formato de entrada

Todos los scripts leen archivos CSV del directorio `combinado/ajustados/` con el formato del pipeline [ddsinaica](https://github.com/JoseAgustin/ddsinaica):

```
eval_<CONT>_<Ciudad>_YYYY-MM-DD.csv
Fecha,Ciudad,max_obs,mod_dia1,mod_dia2,mod_dia3
2026-06-17,Pachuca,45.0617,11.7382,11.3944,13.6731
```

---

## Scripts Python

### 1. `taylor_mensual.py`

**Propósito:** Genera diagramas de Taylor mensuales normalizados para evaluar el desempeño estadístico continuo del modelo (σ, R, CRMSE).

**Metodología:** Cada punto representa Ciudad × Contaminante × Horizonte. Los valores se normalizan por σ_obs para comparar contaminantes con distintas unidades en el mismo diagrama. Implementa la metodología de Taylor (2001, JGR-Atmospheres).

**Control de calidad:** Elimina valores negativos, aplica límites físicos por contaminante y filtro IQR (k=3.0) sobre observaciones de PM2.5.

```bash
# Todos los meses, todas las ciudades
python3 taylor_mensual.py --entrada combinado/ajustados --salida resultados_taylor

# Un mes específico con filtro de ciudades
python3 taylor_mensual.py --mes 2026-06 --ciudades Pachuca Tula CDMX

# Ajustar umbral PM2.5 y factor IQR
python3 taylor_mensual.py --umbral-pm25 120 --iqr-factor 3.0
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--entrada, -i` | Directorio con CSV | `combinado/ajustados` |
| `--salida, -o` | Directorio de salida | `.` |
| `--mes, -m` | Mes específico (YYYY-MM) | todos |
| `--ciudades, -c` | Ciudades a incluir | todas |
| `--umbral-pm25` | Techo absoluto PM2.5 µg/m³ | `500` |
| `--iqr-factor` | Factor k para filtro IQR PM2.5 | `3.0` |
| `--min-pares, -p` | Mínimo pares válidos por serie | `5` |
| `--max-radio, -r` | Radio máximo del diagrama | `1.65` |
| `--dpi` | Resolución de los PNG | `150` |

**Salidas:**
- `taylor_YYYY_MM[_ciudades].png` — diagrama de Taylor
- `estadisticas_taylor[_ciudades].csv` — n, n_orig, n_nan, n_lim, n_iqr, σ_obs, σ_mod, R, BIAS, RMSE, MAE, CRMSE, CRMSE_n, p_valor

---

### 2. `heatmap_desempeno.py`

**Propósito:** Genera mapas de calor (heatmaps) del desempeño temporal del pronóstico con ejes ciudad × mes (o semana), mostrando métricas dicotómicas o continuas con semáforo de colores y bandas de temporada climática.

**Métricas disponibles:** POD, FAR, CSI, TSS, PC, BIAS (dicotómicas) y RMSE, R, BIAS_cont (continuas).

```bash
# POD mensual, todos los contaminantes, umbral categoría Mala
python3 heatmap_desempeno.py --entrada combinado/ajustados --salida resultados_heatmap

# RMSE semanal, solo O3
python3 heatmap_desempeno.py --metrica RMSE --resolucion semana --cont O3

# Categoría Muy Mala, horizonte +48h
python3 heatmap_desempeno.py --metrica POD --categoria muy_mala --horizonte 48h
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--metrica, -M` | POD, FAR, CSI, TSS, PC, BIAS, RMSE, R, BIAS_cont | `POD` |
| `--horizonte, -H` | 24h, 48h, 72h o todos | `todos` |
| `--resolucion, -r` | mes o semana | `mes` |
| `--categoria, -k` | mala o muy_mala (NOM-172) | `mala` |
| `--cont, -C` | O3, PM10, PM25, SO2 | todos |
| `--ciudades, -c` | Ciudades a incluir | todas |
| `--sin-anotar` | Omitir valores numéricos en celdas | — |

**Salidas:**
- `heatmap_<CONT>_<MET>_<HOR>_<RES>_<CAT>.png`

---

### 3. `roebber_desempeno.py`

**Propósito:** Genera el Diagrama de Rendimiento de Roebber (2009) — POD vs SR (= 1−FAR) con isolíneas de CSI (curvas verdes) y BIAS de frecuencia (líneas diagonales). La estrella dorada en (1,1) representa el pronóstico perfecto.

**Modos de agrupación:**
- `horizonte` — colores por ciudad, marcadores por horizonte
- `ciudad` — colores por horizonte
- `mes` — evolución temporal (un punto por mes)
- `temporada` — secas frías / secas calientes / lluvias
- `contaminante` — comparación entre especies

```bash
# Agrupado por ciudad, todos los horizontes
python3 roebber_desempeno.py --cont O3 --agrupar ciudad

# Evolución temporal mes a mes
python3 roebber_desempeno.py --cont O3 --agrupar mes

# Por temporada climática, categoría muy mala
python3 roebber_desempeno.py --agrupar temporada --categoria muy_mala
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--agrupar, -g` | horizonte, ciudad, mes, temporada, contaminante | `horizonte` |
| `--horizonte, -H` | 24h, 48h, 72h o todos | `todos` |
| `--categoria, -k` | mala o muy_mala | `mala` |
| `--separar-cont` | Un PNG independiente por contaminante | — |

**Salidas:**
- `roebber_<CONT>_<AGRUP>_<HOR>_<CAT>.png`

---

### 4. `informe_dicotomico.py`

**Propósito:** Genera un documento Word (.docx) en orientación horizontal (A4 landscape) con estadísticos de verificación dicotómica mensuales (POD, FAR, CSI, TSS, PC, BIAS de frecuencia) y tablas de contingencia (H, M, F, C) por ciudad, contaminante y horizonte. Sin dependencia de Node.js — 100% Python con `python-docx`.

**Umbrales normativos (NOM-172-SEMARNAT-2023):**

| Categoría | O₃ | PM10 | PM2.5 | SO₂ |
|-----------|-----|------|-------|-----|
| **Mala** (naranja) | 135 ppb | 132 µg/m³ | 79 µg/m³ | 185 ppb |
| **Muy Mala** (rojo) | 175 ppb | 213 µg/m³ | 130 µg/m³ | 304 ppb |

```bash
# Categoría mala (default)
python3 informe_dicotomico.py --entrada combinado/ajustados --salida informe_mala.docx

# Categoría muy mala
python3 informe_dicotomico.py --categoria muy_mala --salida informe_muy_mala.docx

# Con CSV de auditoría
python3 informe_dicotomico.py --csv-auditoria stats_dicotomicos.csv
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--categoria, -k` | mala o muy_mala (NOM-172) | `mala` |
| `--mes, -m` | Solo este mes (YYYY-MM) | todos |
| `--ciudades, -c` | Ciudades a incluir | todas |
| `--csv-auditoria` | Exportar CSV con estadísticos | — |
| `--umbral-pm25` | Techo PM2.5 µg/m³ (QC) | `500` |
| `--iqr-factor` | Factor IQR outliers PM2.5 | `3.0` |

**Contenido del documento:**
1. Portada con categoría en color (naranja/rojo)
2. Tabla de umbrales NOM-172 (ambas categorías)
3. Tabla de contingencia conceptual coloreada
4. Definición de métricas con fórmulas
5. Leyenda de semáforo de desempeño
6. Por cada mes → contaminante: tabla de estadísticos + tabla de contingencia H/M/F/C

**Salidas:**
- `informe_dicotomico_<CAT>.docx`
- `dicotomico_stats.csv` (opcional) — N, n_orig, n_nan, n_lim, n_iqr, H, M, F, C, POD, FAR, CSI, TSS, PC, BIAS

---

### 5. `informe_graficas.py`

**Propósito:** Recolecta automáticamente todas las imágenes PNG generadas por `taylor_mensual.py`, `heatmap_desempeno.py` y `roebber_desempeno.py`, genera descripciones técnicas automáticas para cada una, y produce un documento Word (.docx) integrado organizado en tres capítulos.

```bash
# Recolectar de los directorios por defecto
python3 informe_graficas.py

# Especificar directorios y contaminantes
python3 informe_graficas.py \
    --taylor    resultados_taylor \
    --heatmap   resultados_heatmap \
    --roebber   resultados_roebber \
    --salida    informe_visualizaciones.docx \
    --cont      O3 PM10
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--taylor` | Directorio PNG de Taylor | `resultados_taylor` |
| `--heatmap` | Directorio PNG de heatmaps | `resultados_heatmap` |
| `--roebber` | Directorio PNG de Roebber | `resultados_roebber` |
| `--salida, -o` | Ruta del .docx de salida | `informe_visualizaciones.docx` |
| `--cont, -C` | Filtrar por contaminante | todos |

**Contenido del documento:**
- **Cap. 1 — Diagramas de Taylor:** descripción metodológica + imagen + ficha de metadatos por figura
- **Cap. 2 — Mapas de calor:** métrica evaluada, horizonte, umbral y descripción de cómo interpretar cada celda
- **Cap. 3 — Diagramas de Roebber:** agrupación, interpretación de la trayectoria, ubicación de umbrales NOM-172
- **Notas metodológicas:** control de calidad, umbrales, temporadas, horizontes, referencias

---

### 6. `percentiles_obs.py`

**Propósito:** Calcula y visualiza los percentiles empíricos (P10, P25, P50, P75, P90, P95, P99) de las observaciones SINAICA y los compara contra los percentiles del modelo WRF-Chem, organizados por ciudad, mes y contaminante.

**Motivación científica:** Los umbrales normativos nacionales no reflejan la climatología local. El P90 de O₃ en Tula puede diferir en 40+ ppb del P90 de Cuernavaca. Este módulo permite:
1. Conocer la distribución empírica de concentraciones por ciudad y mes
2. Detectar si el modelo subestima sistemáticamente los percentiles extremos
3. Identificar diferencias estacionales en la distribución

```bash
# Todos los contaminantes, P90 de referencia
python3 percentiles_obs.py --entrada combinado/ajustados

# P95, solo O3, ciudades específicas para la serie temporal
python3 percentiles_obs.py --cont O3 --percentil 95 \
    --serie-ciudades CDMX Tula Pachuca

# Con umbral PM2.5 ajustado
python3 percentiles_obs.py --umbral-pm25 120
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--percentil, -p` | P10, P25, P50, P75, P90, P95, P99 | `90` |
| `--serie-ciudades` | Ciudades para serie temporal | primeras 3 |
| `--horizonte, -H` | 24h, 48h, 72h o todos | `todos` |

**Salidas por contaminante:**
- `percentiles_<CONT>_P<N>_heatmap_obs_<HOR>.png` — P_N observado (grande) y modelado (pequeño) por ciudad × mes
- `percentiles_<CONT>_P<N>_heatmap_sesgo_<HOR>.png` — sesgo = (mod/obs − 1), divergente rojo/azul
- `percentiles_<CONT>_QQ_<HOR>.png` — diagrama Q-Q por temporada con marcadores en P50, P90 y P95
- `percentiles_<CONT>_serie_<Ciudad>_<HOR>.png` — serie mensual P25/P50/P75/P90 con bandas IQR
- `tabla_percentiles_<CONT>.csv` — ciudad × mes × horizonte con P10..P99 y métricas derivadas

---

### 7. `curva_habilidad.py`

**Propósito:** Analiza la **sensibilidad del pronóstico WRF-Chem a distintos umbrales de concentración** mediante el barrido continuo de niveles de alerta. Genera la "trayectoria de habilidad" en el espacio de Roebber (2009) al variar el umbral, revelando el comportamiento estructural del modelo independientemente del umbral normativo.

**Preguntas que responde:**
1. ¿En qué umbral el modelo maximiza su habilidad (CSI_max)?
2. ¿El umbral normativo NOM-172 coincide con el óptimo para WRF-Chem?
3. ¿El modelo tiene sesgo de frecuencia independiente del umbral?
4. ¿La habilidad difiere entre temporadas (secas vs lluvias)?

```bash
# Agrupado por ciudad, barrido completo
python3 curva_habilidad.py --cont O3 --agrupar ciudad

# Por temporada climática
python3 curva_habilidad.py --cont O3 --agrupar temporada --horizonte 24h

# Por horizonte (los tres en un diagrama)
python3 curva_habilidad.py --cont O3 --agrupar horizonte

# Paso más fino para publicación
python3 curva_habilidad.py --cont O3 --paso 2 --min-eventos 5
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--agrupar, -g` | ciudad, temporada, horizonte | `ciudad` |
| `--paso` | Paso del barrido de umbrales | `5.0` |
| `--min-eventos` | Mínimo de eventos obs por punto | `3` |
| `--sin-suavizar` | Desactivar suavizado Savitzky-Golay | — |

**Salidas:**
- `curva_habilidad_<CONT>_<HOR>_<AGRUP>.png` — panel 4 subgráficas:
  - (a) POD y SR vs umbral con banda de incertidumbre binomial 90%
  - (b) CSI y TSS vs umbral con marcador del CSI_max
  - (c) BIAS de frecuencia vs umbral (línea ideal en b=1)
  - (d) Trayectoria de habilidad en el espacio de Roebber, coloreada por nivel de umbral; diamantes ◆ = umbrales NOM-172
- `tabla_barrido_umbrales.csv` — umbral × agrupación × horizonte con POD, SR, FAR, CSI, TSS, BIAS, H, M, F, C

---

### 8. `analisis_espacial_percentil.py`

**Propósito:** Realiza la evaluación dicotómica usando **umbrales adaptativos locales** derivados de los percentiles empíricos de las observaciones de cada ciudad, en lugar del umbral normativo nacional fijo. Compara el diagnóstico adaptativo vs el diagnóstico con umbral NOM-172.

**Motivación científica:** Una ciudad con concentraciones basales altas (Tula, fuentes industriales) tendrá muchos eventos incluso con un modelo mediocre si se usa el umbral NOM-172 nacional. El umbral adaptativo por percentil local responde la pregunta: ¿el modelo detecta los días de concentración extrema **para cada ciudad**?

**Tres modos de umbral adaptativo:**

| Modo | Umbral de evento | Pregunta |
|------|-----------------|----------|
| `ciudad_mes` | P_N de obs de esa ciudad en ese mes | ¿Detecta los extremos locales mensuales? |
| `ciudad` | P_N de obs de esa ciudad en todo el período | ¿Detecta los extremos de cada ciudad? |
| `dominio` | P_N de obs de todas las ciudades | ¿Detecta los extremos del dominio? |

```bash
# P90 local por ciudad × mes (más local)
python3 analisis_espacial_percentil.py --cont O3 --percentiles 90

# P75 y P90, modo ciudad, con comparación vs NOM-172
python3 analisis_espacial_percentil.py --cont O3 --percentiles 75 90 \
    --modo ciudad --comparar-nom172

# Modo dominio
python3 analisis_espacial_percentil.py --modo dominio --cont PM10
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--percentiles, -p` | Uno o varios percentiles | `90` |
| `--modo, -m` | ciudad_mes, ciudad, dominio | `ciudad_mes` |
| `--comparar-nom172` | Generar panel comparativo | — |

**Salidas por contaminante × horizonte × percentil:**
- `espacial_<CONT>_<HOR>_<MODO>_P<N>.png` — heatmap 4 métricas (POD/FAR/CSI/TSS) con umbral adaptativo
- `espacial_<CONT>_<HOR>_umbral_P<N>_<MODO>.png` — valor del umbral adaptativo (izq.) y diferencia vs NOM-172 (der.)
- `espacial_<CONT>_<HOR>_delta_P<N>_<MODO>.png` — diferencia métrica_adapt − métrica_nom172 (verde=mejora)
- `espacial_<CONT>_<HOR>_comparacion_P<N>_<MODO>.png` — panel lado a lado adaptativo vs NOM-172
- `tabla_espacial_percentil_<CONT>.csv` — umbral_adapt, umbral_nom172, dif_umbral, métricas _adapt y _nom172, delta_POD/FAR/CSI/TSS

---

### 9. `pronostico_ponderado.py`

**Propósito:** Combina los tres horizontes de pronóstico WRF-Chem (+24h, +48h, +72h) mediante pesos derivados del desempeño histórico (RMSE inverso) para generar un **pronóstico ponderado óptimo para HOY y para MAÑANA**, separado por temporada climática.

**Esquema de horizontes:**
```
Archivo eval_<CONT>_<Ciudad>_FECHA.csv:
  mod_dia1 → run de (FECHA−1) → +24h → válido para FECHA
  mod_dia2 → run de (FECHA−2) → +48h → válido para FECHA
  mod_dia3 → run de (FECHA−3) → +72h → válido para FECHA

HOY   (T):   P = w1·dia1(T) + w2·dia2(T) + w3·dia3(T)
MAÑANA(T+1): P = w1'·dia1(T+1) + w2'·dia2(T)  [pesos renormalizados]
```

**Métodos de ponderación:**

| Método | Descripción |
|--------|-------------|
| `rmse_inv` | w_i = (1/RMSE_i) / Σ(1/RMSE_j) — mayor peso al horizonte más preciso |
| `igual` | w_i = 1/3 — sin preferencia entre horizontes |
| `manual` | `--w1 --w2 --w3` — pesos definidos por el usuario (se renormalizan) |

Los pesos se calculan **por contaminante y por temporada climática** para capturar la variabilidad estacional del desempeño.

```bash
# Pesos por RMSE inverso, fecha automática (último día disponible)
python3 pronostico_ponderado.py --entrada combinado/ajustados

# Fecha específica
python3 pronostico_ponderado.py --fecha 2026-06-30

# Pesos iguales (referencia)
python3 pronostico_ponderado.py --pesos igual

# Pesos manuales
python3 pronostico_ponderado.py --pesos manual --w1 0.5 --w2 0.3 --w3 0.2

# Sin separación por temporada
python3 pronostico_ponderado.py --sin-temporada

# Solo O3 en ciudades específicas
python3 pronostico_ponderado.py --cont O3 --ciudades CDMX Tula Pachuca
```

**Parámetros clave:**

| Parámetro | Descripción | Default |
|-----------|-------------|---------|
| `--fecha, -f` | Fecha de referencia HOY (YYYY-MM-DD) | último día disponible |
| `--pesos, -p` | rmse_inv, igual, manual | `rmse_inv` |
| `--w1, --w2, --w3` | Pesos manuales | 0.5, 0.3, 0.2 |
| `--sin-temporada` | Pesos globales sin separar por temporada | — |
| `--cont, -C` | O3, PM10, PM25, SO2 | todos |
| `--ciudades, -c` | Ciudades a incluir | todas |

**Salidas:**
- `pesos_historicos.csv` — RMSE, MAE, R, n, w1, w2, w3 por contaminante × temporada
- `pronostico_ponderado_HOY.csv` — mod_ponderado_hoy, incertidumbre_1σ, obs_hoy, w1, w2, w3 por ciudad × contaminante
- `pronostico_ponderado_MANANA.csv` — mod_ponderado_manana, horizontes_usados, nota de disponibilidad
- `serie_ponderada_historica.csv` — pronóstico ponderado para todo el histórico (validación)
- `pesos_por_temporada.png` — barras apiladas de pesos con RMSE en eje secundario
- `pronostico_ponderado_<CONT>_serie.png` — serie temporal con banda ±1σ, coloreada por temporada
- `pronostico_ponderado_scatter.png` — scatter obs vs ponderado con regresión y R² por contaminante

---

## Ciudades y contaminantes

### Dominio WRF-Chem

| Ciudad | Red SINAICA | Lat S | Lat N | Lon O | Lon E |
|--------|-------------|-------|-------|-------|-------|
| CDMX | Valle de México | 19.20 | 19.70 | −99.30 | −98.85 |
| Cuernavaca | Cuernavaca | 18.89 | 18.98 | −99.26 | −99.14 |
| Pachuca | Pachuca / Mineral de la Reforma | 20.03 | 20.13 | −98.80 | −98.67 |
| Puebla | Puebla | 18.95 | 19.12 | −98.32 | −98.10 |
| SJdelRio | San Juan del Río | 20.36 | 20.41 | −100.01 | −99.93 |
| Tlaxcala | Tlaxcala | 19.29 | 19.36 | −98.26 | −98.15 |
| Toluca | Toluca | 19.23 | 19.39 | −99.72 | −99.50 |
| Tula | Tula / Tepeji / Atitalaquia | 19.89 | 20.18 | −99.44 | −99.09 |

### Contaminantes evaluados

| Clave | Nombre | Unidad | O₃ | PM10 | PM2.5 | SO₂ |
|-------|--------|--------|----|------|-------|-----|
| `O3` | Ozono | ppbv | ✓ | — | — | — |
| `PM10` | PM10 | µg/m³ | — | ✓ | — | — |
| `PM25` | PM2.5 | µg/m³ | — | — | ✓ | — |
| `SO2` | Dióxido de azufre | ppbv | — | — | — | ✓ (solo Tula) |

---

## Umbrales normativos

**NOM-172-SEMARNAT-2023** (DOF 25/01/2024) — Índice AIRE Y SALUD:

| Categoría | Color | O₃ | PM10 | PM2.5 | SO₂ | Métrica base |
|-----------|-------|-----|------|-------|-----|--------------|
| **Mala** | 🟠 Naranja | 135 ppb | 132 µg/m³ | 79 µg/m³ | 185 ppb | Prom. 1h / prom. móvil pond. 12h |
| **Muy Mala** | 🔴 Rojo | 175 ppb | 213 µg/m³ | 130 µg/m³ | 304 ppb | Ídem |

---

## Temporadas climáticas

Definidas para el centro de México y usadas en todos los módulos:

| Temporada | Meses | Color en gráficas | Característica |
|-----------|-------|-------------------|----------------|
| **Secas frías** | nov–feb | Azul | Menor convección, máx. PM10 |
| **Secas calientes** | mar–may | Naranja | Máximo anual de O₃ fotoquímico |
| **Lluvias** | jun–oct | Verde | Mayor convección, lavado de contaminantes |

---

## Flujo de trabajo recomendado

```bash
# 1. Diagramas de Taylor mensuales (desempeño continuo)
python3 taylor_mensual.py --umbral-pm25 120 --iqr-factor 3.0

# 2. Heatmaps de desempeño (visión temporal)
python3 heatmap_desempeno.py --metrica POD --categoria mala
python3 heatmap_desempeno.py --metrica RMSE

# 3. Diagramas de Roebber (desempeño dicotómico)
python3 roebber_desempeno.py --agrupar temporada
python3 roebber_desempeno.py --agrupar horizonte

# 4. Informe dicotómico Word (ambas categorías NOM-172)
python3 informe_dicotomico.py --categoria mala     --salida informe_mala.docx
python3 informe_dicotomico.py --categoria muy_mala --salida informe_muy_mala.docx

# 5. Informe integrado con todas las visualizaciones
python3 informe_graficas.py --salida informe_visualizaciones.docx

# 6. Análisis de percentiles locales
python3 percentiles_obs.py --percentil 90 --umbral-pm25 120

# 7. Curva de habilidad multi-umbral
python3 curva_habilidad.py --agrupar temporada --cont O3
python3 curva_habilidad.py --agrupar horizonte --cont PM10

# 8. Análisis espacial con umbral adaptativo
python3 analisis_espacial_percentil.py --percentiles 75 90 \
    --modo ciudad_mes --comparar-nom172

# 9. Pronóstico ponderado para hoy y mañana
python3 pronostico_ponderado.py --pesos rmse_inv --umbral-pm25 120
```

---

## Changelog

| Versión | Cambios |
|---------|---------|
| **v1.0.0** | Suite inicial con 9 módulos: taylor_mensual, heatmap_desempeno, roebber_desempeno, informe_dicotomico, informe_graficas, percentiles_obs, curva_habilidad, analisis_espacial_percentil, pronostico_ponderado |

---

## Referencias

- Taylor, K.E., 2001: Summarizing multiple aspects of model performance in a single diagram. *JGR-Atmospheres*, 106(D7), 7183–7192.
- Roebber, P.J., 2009: Visualizing multiple measures of forecast quality. *Wea. Forecasting*, 24, 749–755. doi:10.1175/2008WAF2222159.1
- NOM-172-SEMARNAT-2023: Índice AIRE Y SALUD. DOF 25/01/2024.
- NOM-025-SSA1-2021: Valores límite permisibles para la concentración de PM10 y PM2.5.
- NOM-020-SSA1-2021: Valor límite permisible para la concentración de ozono.
- NOM-022-SSA1-2019: Valor límite permisible para SO₂.

---

## Licencia

MIT License — Copyright © 2026 ICAyCC, UNAM

Código fuente del pipeline operativo: [ddsinaica](https://github.com/JoseAgustin/ddsinaica)
