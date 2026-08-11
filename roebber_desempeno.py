#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
roebber_desempeno.py
====================
Genera el Diagrama de Rendimiento (Performance Diagram) de Roebber (2009)
para la evaluación del pronóstico WRF-Chem contra observaciones SINAICA/INECC.

Fundamento matemático (Roebber, 2009, Wea. Forecasting, 24:749–755)
--------------------------------------------------------------------
A partir de la tabla de contingencia 2×2:

    H  — acierto          (obs EVENTO,    mod EVENTO)
    M  — fallo            (obs EVENTO,    mod NO-evento)
    F  — falsa alarma     (obs NO-evento, mod EVENTO)
    C  — rechazo correcto (obs NO-evento, mod NO-evento)

Se definen:
    POD  = H / (H + M)          Probabilidad de detección       [0, 1]
    SR   = H / (H + F)          Razón de éxito (Success Ratio) [0, 1]
         = 1 − FAR

Relaciones geométricas explotadas en el diagrama:
    CSI  = 1 / (1/POD + 1/SR − 1)     → isolíneas CURVAS
    BIAS = POD / SR                    → isolíneas RECTAS desde el origen

El punto de pronóstico perfecto se ubica en (SR=1, POD=1).
La dirección de mejora óptima (sin cambio de BIAS) es la diagonal a 45°.

Organización de los puntos graficados
--------------------------------------
Cada punto representa una combinación única de:
    Ciudad × Contaminante × Horizonte × [Período opcional]

Opciones de agrupación (--agrupar):
    horizonte   : un punto por ciudad, coloreado por horizonte (+24h/+48h/+72h)
    ciudad      : un punto por horizonte, coloreado por ciudad
    mes         : un punto por mes, coloreado por mes (evolución temporal)
    temporada   : un punto por temporada climática (secas/lluvias)
    contaminante: un punto por contaminante

Umbrales normativos — NOM-172-SEMARNAT-2023
--------------------------------------------
    Categoría "mala"     : O3=135ppb  PM10=132µg/m³  PM25=79µg/m³  SO2=185ppb
    Categoría "muy_mala" : O3=175ppb  PM10=213µg/m³  PM25=130µg/m³ SO2=304ppb

Uso
---
    # Diagrama básico: todos los contaminantes, agrupar por horizonte
    python3 roebber_desempeno.py

    # Solo O3, agrupar por ciudad
    python3 roebber_desempeno.py --cont O3 --agrupar ciudad

    # Evolución temporal: un punto por mes
    python3 roebber_desempeno.py --cont PM10 --agrupar mes --ciudades Tula

    # Comparar temporadas
    python3 roebber_desempeno.py --agrupar temporada --cont O3

    # Categoría muy mala, horizonte +48h
    python3 roebber_desempeno.py --categoria muy_mala --horizonte 48h

    # Un diagrama por contaminante
    python3 roebber_desempeno.py --separar-cont

    python3 roebber_desempeno.py --help

Salidas
-------
    roebber_<CONT>_<AGRUPAR>_<HOR>_<CATEGORIA>.png  — un PNG por diagrama

Dependencias
------------
    pip install pandas numpy matplotlib scipy

Referencias
-----------
    Roebber, P. J., 2009: Visualizing multiple measures of forecast quality.
        Wea. Forecasting, 24, 749–755. doi:10.1175/2008WAF2222159.1

Autor  : Pipeline ddsinaica / WRF-Chem — ICAyCC, UNAM
Versión: 1.0.0 (2026-07)
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from scipy import stats

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────

CATEGORIAS_NOM172 = {
    "mala": {
        "O3":   135.0,   # ppb
        "PM10": 132.0,   # µg/m³
        "PM25":  79.0,   # µg/m³
        "SO2":  185.0,   # ppb
    },
    "muy_mala": {
        "O3":   175.0,
        "PM10": 213.0,
        "PM25": 130.0,
        "SO2":  304.0,
    },
}
ALIAS_CATEGORIA = {
    "mala": "mala", "malo": "mala", "alto": "mala", "naranja": "mala",
    "muy_mala": "muy_mala", "muymala": "muy_mala",
    "muy_alto": "muy_mala", "rojo": "muy_mala",
}

CIUDADES_DOMINIO = [
    "CDMX", "Cuernavaca", "Pachuca", "Puebla",
    "SJdelRio", "Tlaxcala", "Toluca", "Tula",
]
_CIUDADES_NORM = {c.lower(): c for c in CIUDADES_DOMINIO}
_CIUDADES_NORM.update({
    "sjdelrio": "SJdelRio", "san juan del rio": "SJdelRio",
    "cdmx": "CDMX", "ciudad de mexico": "CDMX",
})

CONTAMINANTES_ORDEN = ["O3", "PM10", "PM25", "SO2"]
META_CONT = {
    "O3":   {"nombre": "Ozono (O₃)",               "unidad": "ppbv"},
    "PM10": {"nombre": "PM10",                      "unidad": "µg/m³"},
    "PM25": {"nombre": "PM2.5",                     "unidad": "µg/m³"},
    "SO2":  {"nombre": "Dióxido de azufre (SO₂)",  "unidad": "ppbv"},
}

HORIZONTES_MAP = {
    "24h": "mod_dia1", "+24h": "mod_dia1", "dia1": "mod_dia1",
    "48h": "mod_dia2", "+48h": "mod_dia2", "dia2": "mod_dia2",
    "72h": "mod_dia3", "+72h": "mod_dia3", "dia3": "mod_dia3",
    "todos": None,
}
HORIZONTES_COLS   = ["mod_dia1", "mod_dia2", "mod_dia3"]
HORIZONTES_LABELS = {
    "mod_dia1": "+24 h", "mod_dia2": "+48 h", "mod_dia3": "+72 h",
}

# Control de calidad (mismo criterio que el resto del pipeline)
LIMITES_VALIDOS = {
    "O3":   (0.0,  300.0,  0.0,  300.0),
    "PM10": (0.0, 1000.0,  0.0, 1000.0),
    "PM25": (0.0,  500.0,  0.0,  500.0),
    "SO2":  (0.0, 1000.0,  0.0, 1000.0),
}
CONTAMINANTES_FILTRO_IQR = {"PM25"}
IQR_FACTOR  = 3.0
UMBRAL_PM25 = 500.0
MIN_DIAS    = 5

# Temporadas climáticas del centro de México
TEMPORADAS = {
    "Secas frías":      [11, 12, 1, 2],
    "Secas calientes":  [3, 4, 5],
    "Lluvias":          [6, 7, 8, 9, 10],
}
COLOR_TEMPORADA = {
    "Secas frías":      "#4A90D9",
    "Secas calientes":  "#F5A623",
    "Lluvias":          "#7ED321",
}

NOMBRE_MES = {
    1: "Ene", 2: "Feb",  3: "Mar",  4: "Abr",
    5: "May", 6: "Jun",  7: "Jul",  8: "Ago",
    9: "Sep", 10: "Oct", 11: "Nov", 12: "Dic",
}

# Paletas de colores
COLORES_HOR = {
    "mod_dia1": "#1f77b4",   # azul  +24h
    "mod_dia2": "#ff7f0e",   # naranja +48h
    "mod_dia3": "#d62728",   # rojo  +72h
}
MARKERS_HOR = {
    "mod_dia1": "o",
    "mod_dia2": "s",
    "mod_dia3": "^",
}
COLORES_CONT = {
    "O3":   "#1f77b4",
    "PM10": "#d62728",
    "PM25": "#ff7f0e",
    "SO2":  "#9467bd",
}
# Paleta de 8 colores para ciudades
COLORES_CIUDAD = {
    "CDMX":       "#e6194b",
    "Cuernavaca": "#f58231",
    "Pachuca":    "#ffe119",
    "Puebla":     "#3cb44b",
    "SJdelRio":   "#42d4f4",
    "Tlaxcala":   "#4363d8",
    "Toluca":     "#911eb4",
    "Tula":       "#a9a9a9",
}

DPI = 150


# ──────────────────────────────────────────────────────────────────────────────
# UTILIDADES
# ──────────────────────────────────────────────────────────────────────────────

def normalizar_ciudad(nombre: str) -> Optional[str]:
    return _CIUDADES_NORM.get(nombre.strip().lower())

def parsear_lista_ciudades(valor) -> Optional[List[str]]:
    if not valor:
        return None
    crudos = []
    items = [valor] if isinstance(valor, str) else valor
    for item in items:
        crudos.extend(t for t in re.split(r"[,\s]+", item.strip()) if t)
    canonicas, invalidas = [], []
    for c in crudos:
        cn = normalizar_ciudad(c)
        if cn is None: invalidas.append(c)
        elif cn not in canonicas: canonicas.append(cn)
    if invalidas:
        sys.exit(f"[ERROR] Ciudad(es) no reconocida(s): {invalidas}\n"
                 f"        Válidas: {CIUDADES_DOMINIO}")
    return canonicas

def parsear_nombre(ruta: str) -> Optional[Tuple[str, str]]:
    m = re.match(r"^eval_([A-Za-z0-9]+)_(.+)_\d{4}-\d{2}-\d{2}$", Path(ruta).stem)
    if not m: return None
    cont   = m.group(1).upper()
    ciudad = normalizar_ciudad(m.group(2)) or m.group(2)
    return cont, ciudad

def mes_a_temporada(mes: int) -> str:
    for nombre, meses in TEMPORADAS.items():
        if mes in meses:
            return nombre
    return "Lluvias"


# ──────────────────────────────────────────────────────────────────────────────
# CONTROL DE CALIDAD
# ──────────────────────────────────────────────────────────────────────────────

def limpiar_serie(
    obs: np.ndarray,
    mod: np.ndarray,
    cont: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Aplica QC: NaN, límites físicos y filtro IQR para PM2.5."""
    mask = np.isfinite(obs) & np.isfinite(mod)
    obs, mod = obs[mask], mod[mask]

    lim = LIMITES_VALIDOS.get(cont)
    if lim:
        min_obs, max_obs, min_mod, max_mod = lim
        if cont == "PM25":
            max_obs = UMBRAL_PM25
        mask_lim = (
            (obs >= min_obs) & (obs <= max_obs) &
            (mod >= min_mod) & (mod <= max_mod)
        )
        obs, mod = obs[mask_lim], mod[mask_lim]

    if cont in CONTAMINANTES_FILTRO_IQR and len(obs) >= 4:
        q1, q3 = np.percentile(obs, [25, 75])
        iqr = q3 - q1
        if iqr > 0:
            mask_iqr = (obs >= q1 - IQR_FACTOR*iqr) & (obs <= q3 + IQR_FACTOR*iqr)
            obs, mod = obs[mask_iqr], mod[mask_iqr]

    return obs, mod


# ──────────────────────────────────────────────────────────────────────────────
# CÁLCULO DE POD Y SR
# ──────────────────────────────────────────────────────────────────────────────

def calcular_pod_sr(
    obs: np.ndarray,
    mod: np.ndarray,
    umbral: float,
) -> Optional[Tuple[float, float, int]]:
    """
    Calcula POD y SR (Success Ratio = 1 - FAR) a partir de arrays ya limpios.

    Retorna (POD, SR, N) o None si N < MIN_DIAS.
    """
    if len(obs) < MIN_DIAS:
        return None

    obs_ev = obs >= umbral
    mod_ev = mod >= umbral
    H = int(np.sum( obs_ev &  mod_ev))
    M = int(np.sum( obs_ev & ~mod_ev))
    F = int(np.sum(~obs_ev &  mod_ev))

    def safe(num, den): return num / den if den > 0 else np.nan

    POD = safe(H, H + M)
    SR  = safe(H, H + F)     # = 1 - FAR

    if not (np.isfinite(POD) and np.isfinite(SR)):
        return None

    N = len(obs)
    return float(POD), float(SR), N


# ──────────────────────────────────────────────────────────────────────────────
# LECTURA DE DATOS
# ──────────────────────────────────────────────────────────────────────────────

def leer_datos(
    directorio: str,
    ciudades_filtro: Optional[List[str]],
    inicio: Optional[str],
    fin: Optional[str],
    cont_filtro: Optional[str],
) -> pd.DataFrame:
    """Lee los eval_*.csv y devuelve un DataFrame largo con columnas normalizadas."""
    archivos = sorted(glob.glob(os.path.join(directorio, "eval_*.csv")))
    if not archivos:
        sys.exit(f"[ERROR] No se encontraron eval_*.csv en '{directorio}'.")

    partes = []
    for ruta in archivos:
        meta = parsear_nombre(ruta)
        if not meta: continue
        cont, ciudad = meta
        if cont not in CATEGORIAS_NOM172["mala"]: continue
        if cont_filtro and cont != cont_filtro: continue
        if ciudades_filtro and ciudad not in ciudades_filtro: continue
        try:
            df = pd.read_csv(ruta, parse_dates=["Fecha"])
        except Exception:
            continue
        if not {"Fecha","max_obs","mod_dia1","mod_dia2","mod_dia3"}.issubset(df.columns):
            continue
        df["ciudad"] = ciudad
        df["contaminante"] = cont
        partes.append(df[["Fecha","ciudad","contaminante",
                           "max_obs","mod_dia1","mod_dia2","mod_dia3"]])

    if not partes:
        sys.exit("[ERROR] Ningún archivo coincidió con los filtros.")

    datos = pd.concat(partes, ignore_index=True).rename(columns={"Fecha": "fecha"})
    datos["fecha"] = pd.to_datetime(datos["fecha"])

    if inicio:
        datos = datos[datos["fecha"] >= pd.Timestamp(inicio + "-01")]
    if fin:
        datos = datos[datos["fecha"] <= pd.Timestamp(fin + "-01") + pd.offsets.MonthEnd(0)]

    print(f"[INFO] Registros cargados: {len(datos):,}  "
          f"({datos['fecha'].dt.date.nunique()} días únicos)")
    return datos


# ──────────────────────────────────────────────────────────────────────────────
# CONSTRUCCIÓN DE PUNTOS DEL DIAGRAMA
# ──────────────────────────────────────────────────────────────────────────────

def construir_puntos(
    datos: pd.DataFrame,
    cont: str,
    umbral: float,
    cols_mod: List[str],
    agrupar: str,
) -> List[Dict]:
    """
    Genera la lista de puntos a graficar.

    Cada punto es un dict con:
        pod, sr, n, etiqueta, color, marker, size, alpha, zorder
    """
    df = datos[datos["contaminante"] == cont].copy()
    if df.empty:
        return []

    df["mes"]      = df["fecha"].dt.month
    df["mes_str"]  = df["fecha"].dt.to_period("M").astype(str)
    df["temporada"] = df["mes"].map(mes_a_temporada)

    puntos = []

    # ── Paleta de meses (colormap continuo) ──────────────────────────────────
    meses_unicos = sorted(df["mes_str"].unique())
    cmap_meses   = matplotlib.colormaps.get_cmap("hsv").resampled(len(meses_unicos) + 1)
    color_mes    = {m: cmap_meses(i) for i, m in enumerate(meses_unicos)}

    # ── Colores de ciudades presentes ────────────────────────────────────────
    ciudades_presentes = sorted(df["ciudad"].unique(),
                                key=lambda c: CIUDADES_DOMINIO.index(c)
                                if c in CIUDADES_DOMINIO else 99)

    for col_mod in cols_mod:
        hor_label = HORIZONTES_LABELS[col_mod]

        if agrupar == "horizonte":
            # Un punto por ciudad (todos los meses juntos)
            for ciudad in ciudades_presentes:
                dfc = df[df["ciudad"] == ciudad]
                obs = np.asarray(dfc["max_obs"].values, float)
                mod = np.asarray(dfc[col_mod].values, float)
                obs, mod = limpiar_serie(obs, mod, cont)
                res = calcular_pod_sr(obs, mod, umbral)
                if res is None: continue
                pod, sr, n = res
                puntos.append({
                    "pod":     pod, "sr": sr, "n": n,
                    "etiqueta": f"{ciudad}\n{hor_label}",
                    "color":   COLORES_CIUDAD.get(ciudad, "gray"),
                    "marker":  MARKERS_HOR[col_mod],
                    "size":    90, "alpha": 0.88, "zorder": 5,
                    "ciudad":  ciudad, "horizonte": hor_label, "grupo": ciudad,
                })

        elif agrupar == "ciudad":
            # Un punto por horizonte (todos los meses juntos)
            for ciudad in ciudades_presentes:
                dfc = df[df["ciudad"] == ciudad]
                obs = np.asarray(dfc["max_obs"].values, float)
                mod = np.asarray(dfc[col_mod].values, float)
                obs, mod = limpiar_serie(obs, mod, cont)
                res = calcular_pod_sr(obs, mod, umbral)
                if res is None: continue
                pod, sr, n = res
                puntos.append({
                    "pod": pod, "sr": sr, "n": n,
                    "etiqueta": f"{ciudad}\n{hor_label}",
                    "color":   COLORES_HOR[col_mod],
                    "marker":  MARKERS_HOR[col_mod],
                    "size":    90, "alpha": 0.88, "zorder": 5,
                    "ciudad":  ciudad, "horizonte": hor_label, "grupo": hor_label,
                })

        elif agrupar == "mes":
            # Un punto por mes
            for mes_str in meses_unicos:
                dfm = df[df["mes_str"] == mes_str]
                for ciudad in ciudades_presentes:
                    dfc = dfm[dfm["ciudad"] == ciudad]
                    if dfc.empty: continue
                    obs = np.asarray(dfc["max_obs"].values, float)
                    mod = np.asarray(dfc[col_mod].values, float)
                    obs, mod = limpiar_serie(obs, mod, cont)
                    res = calcular_pod_sr(obs, mod, umbral)
                    if res is None: continue
                    pod, sr, n = res
                    mm = int(mes_str.split("-")[1])
                    puntos.append({
                        "pod": pod, "sr": sr, "n": n,
                        "etiqueta": f"{NOMBRE_MES[mm]}\n{ciudad[:4]}",
                        "color":   color_mes[mes_str],
                        "marker":  MARKERS_HOR[col_mod],
                        "size":    80, "alpha": 0.85, "zorder": 5,
                        "ciudad":  ciudad, "horizonte": hor_label,
                        "grupo":   mes_str, "mes_num": mm,
                    })

        elif agrupar == "temporada":
            # Un punto por temporada × ciudad
            for temporada, col_temp in COLOR_TEMPORADA.items():
                dft = df[df["temporada"] == temporada]
                for ciudad in ciudades_presentes:
                    dfc = dft[dft["ciudad"] == ciudad]
                    if dfc.empty: continue
                    obs = np.asarray(dfc["max_obs"].values, float)
                    mod = np.asarray(dfc[col_mod].values, float)
                    obs, mod = limpiar_serie(obs, mod, cont)
                    res = calcular_pod_sr(obs, mod, umbral)
                    if res is None: continue
                    pod, sr, n = res
                    puntos.append({
                        "pod": pod, "sr": sr, "n": n,
                        "etiqueta": f"{temporada.split()[0]}\n{ciudad[:4]}",
                        "color":   col_temp,
                        "marker":  MARKERS_HOR[col_mod],
                        "size":    100, "alpha": 0.88, "zorder": 5,
                        "ciudad":  ciudad, "horizonte": hor_label,
                        "grupo":   temporada,
                    })

        elif agrupar == "contaminante":
            # Un punto por contaminante (llamado desde fuera con cont variable)
            for ciudad in ciudades_presentes:
                dfc = df[df["ciudad"] == ciudad]
                obs = np.asarray(dfc["max_obs"].values, float)
                mod = np.asarray(dfc[col_mod].values, float)
                obs, mod = limpiar_serie(obs, mod, cont)
                res = calcular_pod_sr(obs, mod, umbral)
                if res is None: continue
                pod, sr, n = res
                puntos.append({
                    "pod": pod, "sr": sr, "n": n,
                    "etiqueta": f"{cont}\n{ciudad[:4]}",
                    "color":   COLORES_CONT.get(cont, "gray"),
                    "marker":  MARKERS_HOR[col_mod],
                    "size":    90, "alpha": 0.88, "zorder": 5,
                    "ciudad":  ciudad, "horizonte": hor_label,
                    "grupo":   cont,
                })

    return puntos


# ──────────────────────────────────────────────────────────────────────────────
# FONDO DEL DIAGRAMA DE ROEBBER
# ──────────────────────────────────────────────────────────────────────────────

def dibujar_fondo_roebber(ax: plt.Axes):
    """
    Dibuja sobre `ax` el fondo del diagrama de Roebber (2009):
      - Isolíneas de CSI (curvas, en verde)
      - Líneas de BIAS (diagonales punteadas, en gris)
      - Eje perfecto (SR=1, POD=1) marcado con estrella
      - Relleno de fondo por nivel de CSI (gradiente suave)
    """
    SR_grid  = np.linspace(0.001, 1.0, 500)
    POD_grid = np.linspace(0.001, 1.0, 500)
    SR_2d, POD_2d = np.meshgrid(SR_grid, POD_grid)

    # CSI = 1 / (1/SR + 1/POD - 1)
    CSI_2d = 1.0 / (1.0/SR_2d + 1.0/POD_2d - 1.0)
    CSI_2d = np.clip(CSI_2d, 0, 1)

    # ── Relleno de fondo por CSI (muy suave) ─────────────────────────────────
    ax.contourf(SR_2d, POD_2d, CSI_2d,
                levels=np.linspace(0, 1, 21),
                cmap="YlGn", alpha=0.18, zorder=0)

    # ── Isolíneas de CSI ─────────────────────────────────────────────────────
    niveles_csi = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    cs = ax.contour(SR_2d, POD_2d, CSI_2d,
                    levels=niveles_csi,
                    colors="#2d6a2d", linewidths=0.9, alpha=0.7, zorder=1)
    ax.clabel(cs, fmt="CSI=%.1f", fontsize=7.5, inline=True,
              inline_spacing=4, colors="#2d6a2d")

    # ── Líneas de BIAS = POD / SR ─────────────────────────────────────────────
    bias_vals  = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    bias_estilos = {1.0: {"lw": 1.6, "ls": "--", "color": "#333333"},}   # BIAS=1 más grueso
    sr_line = np.linspace(0, 1, 300)

    for bv in bias_vals:
        pod_line = bv * sr_line
        mask_b   = (pod_line >= 0) & (pod_line <= 1.0)
        lw       = 1.6 if bv == 1.0 else 0.85
        ls       = "--"
        col      = "#333333" if bv == 1.0 else "#888888"
        ax.plot(sr_line[mask_b], pod_line[mask_b],
                ls=ls, lw=lw, color=col, alpha=0.65, zorder=1)

        # Etiqueta al final de la línea visible
        idx_end = np.where(mask_b)[0]
        if len(idx_end) == 0: continue
        ie = idx_end[-1]
        x_lbl = sr_line[ie]
        y_lbl = pod_line[ie]
        if x_lbl > 0.02 and y_lbl <= 1.0:
            bias_txt = f"b={bv:.2g}" if bv != 1.0 else "b=1 (no sesgo)"
            ax.text(x_lbl + 0.01, y_lbl, bias_txt,
                    fontsize=6.5, color=col, va="bottom",
                    ha="left", rotation=0, zorder=2,
                    alpha=0.85)

    # ── Punto perfecto ────────────────────────────────────────────────────────
    ax.plot(1.0, 1.0, "*", ms=16, color="gold",
            markeredgecolor="#333333", markeredgewidth=0.8,
            zorder=8, label="Pronóstico perfecto")
    ax.annotate("Pronóstico\nperfecto", xy=(1.0, 1.0),
                xytext=(0.84, 0.92),
                fontsize=7.5, color="#333333",
                arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8))


# ──────────────────────────────────────────────────────────────────────────────
# GRAFICAR EL DIAGRAMA
# ──────────────────────────────────────────────────────────────────────────────

def graficar_roebber(
    puntos: List[Dict],
    titulo: str,
    agrupar: str,
    dir_salida: str,
    nombre_archivo: str,
    anotar: bool = True,
):
    """Dibuja y guarda un diagrama de Roebber completo."""

    fig, ax = plt.subplots(figsize=(8.5, 8))

    # ── Fondo del diagrama ────────────────────────────────────────────────────
    dibujar_fondo_roebber(ax)

    # ── Puntos de datos ───────────────────────────────────────────────────────
    grupos_vistos: Dict[str, Dict] = {}

    for pt in puntos:
        grp = pt.get("grupo", pt["etiqueta"])
        ax.scatter(
            pt["sr"], pt["pod"],
            c=pt["color"], marker=pt["marker"],
            s=pt["size"], alpha=pt["alpha"],
            edgecolors="white", linewidths=0.7,
            zorder=pt["zorder"],
        )
        # Anotación de etiqueta
        if anotar:
            ax.annotate(
                pt["etiqueta"],
                xy=(pt["sr"], pt["pod"]),
                xytext=(pt["sr"] + 0.012, pt["pod"] + 0.008),
                fontsize=6.5, color=pt["color"], zorder=6,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
            )
        # Acumular para leyenda (solo un handle por grupo)
        if grp not in grupos_vistos:
            grupos_vistos[grp] = pt

    # ── Leyenda ───────────────────────────────────────────────────────────────
    handles = []
    if agrupar == "horizonte":
        for col_mod, lbl in HORIZONTES_LABELS.items():
            h = ax.scatter([], [], c="gray", marker=MARKERS_HOR[col_mod],
                           s=60, label=f"Marcador: {lbl}")
            handles.append(h)
        # Leyenda de colores = ciudades
        for grp, pt in grupos_vistos.items():
            h = ax.scatter([], [], c=pt["color"], marker="o",
                           s=60, label=pt.get("ciudad", grp))
            handles.append(h)

    elif agrupar == "ciudad":
        for col_mod, lbl in HORIZONTES_LABELS.items():
            h = ax.scatter([], [], c=COLORES_HOR[col_mod],
                           marker=MARKERS_HOR[col_mod], s=60, label=lbl)
            handles.append(h)

    elif agrupar == "temporada":
        for temp, col in COLOR_TEMPORADA.items():
            h = ax.scatter([], [], c=col, marker="o", s=60, label=temp)
            handles.append(h)
        for col_mod, lbl in HORIZONTES_LABELS.items():
            h = ax.scatter([], [], c="gray", marker=MARKERS_HOR[col_mod],
                           s=50, label=f"Marcador: {lbl}")
            handles.append(h)

    elif agrupar == "mes":
        # Muestra solo los meses que aparecen
        meses_en_puntos = sorted({pt.get("mes_num", 1) for pt in puntos})
        # Para la leyenda de colores necesitamos reconstruir el mapa
        meses_unicos_str = sorted({
            pt.get("grupo", "") for pt in puntos
            if re.match(r"\d{4}-\d{2}", pt.get("grupo",""))
        })
        cmap_m = matplotlib.colormaps.get_cmap("hsv").resampled(len(meses_unicos_str) + 1)
        for i, ms in enumerate(meses_unicos_str):
            mm = int(ms.split("-")[1])
            h = ax.scatter([], [], c=[cmap_m(i)], marker="o",
                           s=55, label=f"{NOMBRE_MES[mm]} {ms[:4]}")
            handles.append(h)
        for col_mod, lbl in HORIZONTES_LABELS.items():
            h = ax.scatter([], [], c="gray", marker=MARKERS_HOR[col_mod],
                           s=50, label=f"Marcador: {lbl}")
            handles.append(h)

    elif agrupar == "contaminante":
        for cont_k, col in COLORES_CONT.items():
            if any(pt.get("grupo") == cont_k for pt in puntos):
                h = ax.scatter([], [], c=col, marker="o",
                               s=60, label=cont_k)
                handles.append(h)

    # Añadir punto perfecto a la leyenda
    handles.append(plt.Line2D([0], [0], marker="*", color="w",
                               markerfacecolor="gold",
                               markeredgecolor="#333333",
                               markersize=11,
                               label="Pronóstico perfecto"))

    ax.legend(handles=handles, loc="lower left",
              fontsize=7.5, framealpha=0.92,
              edgecolor="gray", ncol=2,
              title="Leyenda", title_fontsize=8)

    # ── Ejes y formato ────────────────────────────────────────────────────────
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(-0.02, 1.05)
    ax.set_aspect("equal")
    ax.set_xlabel("SR — Razón de Éxito  (1 − FAR = H / (H + F))",
                  fontsize=10, labelpad=8)
    ax.set_ylabel("POD — Probabilidad de Detección  (H / (H + M))",
                  fontsize=10, labelpad=8)

    ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.grid(True, which="major", color="#CCCCCC", lw=0.5, zorder=0)
    ax.tick_params(labelsize=8.5)

    # ── Título y referencias ──────────────────────────────────────────────────
    ax.set_title(titulo, fontsize=11, fontweight="bold", pad=12)

    fig.text(
        0.01, 0.005,
        "Roebber (2009) Wea. Forecasting 24:749–755  |  "
        "NOM-172-SEMARNAT-2023  |  WRF-Chem vs SINAICA/INECC  |  ICAyCC, UNAM",
        fontsize=6.5, style="italic", color="gray",
    )

    # Nota de interpretación
    fig.text(
        0.99, 0.005,
        "↗ mejora óptima (45°)  |  líneas — = BIAS  |  curvas — = CSI",
        fontsize=6.5, color="#444444", ha="right",
    )

    plt.tight_layout(rect=[0, 0.02, 1, 1])

    os.makedirs(dir_salida, exist_ok=True)
    ruta = os.path.join(dir_salida, nombre_archivo)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# ARGUMENTOS CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="roebber_desempeno.py",
        description=(
            "Genera el Diagrama de Rendimiento de Roebber (2009): "
            "POD vs SR con isolíneas de CSI y líneas de BIAS, "
            "para la evaluación WRF-Chem vs SINAICA/INECC."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--entrada", "-i", default="combinado/ajustados",
                   help="Directorio con eval_*.csv (default: combinado/ajustados)")
    p.add_argument("--salida",  "-o", default="resultados_roebber",
                   help="Directorio de salida PNG (default: resultados_roebber)")
    p.add_argument("--cont", "-C", default=None, choices=CONTAMINANTES_ORDEN,
                   help="Contaminante a procesar (default: todos en diagramas separados)")
    p.add_argument("--horizonte", "-H", default="todos",
                   help="Horizonte: 24h, 48h, 72h o 'todos' (default: todos)")
    p.add_argument("--agrupar", "-g", default="horizonte",
                   choices=["horizonte", "ciudad", "mes", "temporada", "contaminante"],
                   help=(
                       "Cómo colorear/agrupar los puntos:\n"
                       "  horizonte   — colores=ciudad, marcadores=horizonte (default)\n"
                       "  ciudad      — colores=horizonte, un punto por ciudad\n"
                       "  mes         — un punto por mes (evolución temporal)\n"
                       "  temporada   — un punto por temporada climática\n"
                       "  contaminante— colores=contaminante (un diagrama conjunto)"
                   ))
    p.add_argument("--separar-cont", action="store_true",
                   help="Generar un PNG independiente por cada contaminante")
    p.add_argument("--categoria", "-k", default="mala",
                   help="Umbral NOM-172: 'mala' o 'muy_mala' (default: mala)")
    p.add_argument("--ciudades", "-c", nargs="+", default=None,
                   help=f"Ciudades a incluir. Catálogo: {', '.join(CIUDADES_DOMINIO)}")
    p.add_argument("--inicio", default=None, metavar="YYYY-MM",
                   help="Mes de inicio del período (ej: 2026-01)")
    p.add_argument("--fin",    default=None, metavar="YYYY-MM",
                   help="Mes de fin del período (ej: 2026-12)")
    p.add_argument("--umbral-pm25", type=float, default=UMBRAL_PM25,
                   help=f"Techo absoluto PM2.5 µg/m³ (default: {UMBRAL_PM25})")
    p.add_argument("--iqr-factor",  type=float, default=IQR_FACTOR,
                   help=f"Factor IQR outliers PM2.5 (default: {IQR_FACTOR})")
    p.add_argument("--min-dias",    type=int,   default=MIN_DIAS,
                   help=f"Mínimo días válidos por punto (default: {MIN_DIAS})")
    p.add_argument("--sin-anotar", action="store_true",
                   help="Omitir etiquetas sobre los puntos")
    p.add_argument("--dpi", type=int, default=DPI,
                   help=f"Resolución de los PNG (default: {DPI})")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# PUNTO DE ENTRADA
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    global UMBRAL_PM25, IQR_FACTOR, MIN_DIAS, DPI
    UMBRAL_PM25 = args.umbral_pm25
    IQR_FACTOR  = args.iqr_factor
    MIN_DIAS    = args.min_dias
    DPI         = args.dpi

    # Resolver categoría
    cat_key = ALIAS_CATEGORIA.get(args.categoria.lower().replace(" ", "_"))
    if cat_key is None:
        sys.exit(f"[ERROR] Categoría no reconocida: '{args.categoria}'\n"
                 f"        Opciones: mala, muy_mala")
    cat_label = "Mala (naranja)" if cat_key == "mala" else "Muy Mala (rojo)"

    # Resolver horizontes
    if args.horizonte.lower() in ("todos", "all"):
        cols_mod = HORIZONTES_COLS
        hor_txt  = "todos los horizontes"
        hor_slug = "todos"
    else:
        col = HORIZONTES_MAP.get(args.horizonte.lower().replace("+", ""))
        if col is None:
            sys.exit(f"[ERROR] Horizonte no reconocido: '{args.horizonte}'")
        cols_mod = [col]
        hor_txt  = HORIZONTES_LABELS[col]
        hor_slug = hor_txt.replace(" ", "")

    ciudades_filtro = parsear_lista_ciudades(args.ciudades)

    # Contaminantes a procesar
    conts = [args.cont] if args.cont else CONTAMINANTES_ORDEN

    print("=" * 66)
    print("  Diagrama de Roebber — WRF-Chem / SINAICA")
    print("=" * 66)
    print(f"  Entrada      : {args.entrada}")
    print(f"  Salida       : {args.salida}")
    print(f"  Contaminante : {conts}")
    print(f"  Horizonte    : {hor_txt}")
    print(f"  Agrupación   : {args.agrupar}")
    print(f"  Categoría    : {cat_key} ({cat_label})")
    print(f"  Ciudades     : {ciudades_filtro or 'todas'}")
    print(f"  Período      : {args.inicio or 'inicio'} – {args.fin or 'fin'}")
    print("=" * 66)

    # Leer datos una sola vez
    datos = leer_datos(
        args.entrada, ciudades_filtro,
        args.inicio, args.fin,
        args.cont if args.separar_cont or args.cont else None,
    )

    # ── Modo "contaminante" en un solo diagrama ───────────────────────────────
    if args.agrupar == "contaminante" and not args.separar_cont:
        todos_puntos = []
        for cont in conts:
            umbral = CATEGORIAS_NOM172[cat_key][cont]
            todos_puntos += construir_puntos(
                datos, cont, umbral, cols_mod, "contaminante"
            )
        if not todos_puntos:
            sys.exit("[ERROR] Sin puntos válidos para graficar.")

        titulo = (
            f"Diagrama de Roebber — Todos los contaminantes\n"
            f"Horizonte: {hor_txt}  |  Umbral: {cat_label} (NOM-172-SEMARNAT-2023)"
        )
        nombre_png = f"roebber_todos_{args.agrupar}_{hor_slug}_{cat_key}.png"
        graficar_roebber(todos_puntos, titulo, "contaminante",
                         args.salida, nombre_png, not args.sin_anotar)

    # ── Modo separado: un PNG por contaminante ────────────────────────────────
    else:
        for cont in conts:
            if cont not in datos["contaminante"].unique():
                print(f"[WARN] {cont}: sin datos — omitido.")
                continue

            umbral = CATEGORIAS_NOM172[cat_key][cont]
            meta   = META_CONT[cont]
            print(f"\n── {cont} ({cat_key}: {umbral} {meta['unidad']}) ──")

            puntos = construir_puntos(
                datos, cont, umbral, cols_mod, args.agrupar
            )
            if not puntos:
                print(f"[WARN] {cont}: sin puntos válidos (¿excedencias insuficientes?).")
                continue

            n_pts = len(puntos)
            pod_vals = [pt["pod"] for pt in puntos]
            sr_vals  = [pt["sr"]  for pt in puntos]
            print(f"       {n_pts} puntos  "
                  f"POD=[{min(pod_vals):.2f}–{max(pod_vals):.2f}]  "
                  f"SR=[{min(sr_vals):.2f}–{max(sr_vals):.2f}]")

            titulo = (
                f"Diagrama de Roebber — {meta['nombre']}\n"
                f"Horizonte: {hor_txt}  |  "
                f"Umbral {cat_label}: {umbral} {meta['unidad']}  "
                f"(NOM-172-SEMARNAT-2023)"
            )
            nombre_png = (
                f"roebber_{cont}_{args.agrupar}_{hor_slug}_{cat_key}.png"
            )
            graficar_roebber(puntos, titulo, args.agrupar,
                             args.salida, nombre_png, not args.sin_anotar)

    print(f"\n[DONE] Diagramas guardados en '{args.salida}/'")


if __name__ == "__main__":
    main()

