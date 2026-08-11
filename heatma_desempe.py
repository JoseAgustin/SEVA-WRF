#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
heatmap_desempeno.py
====================
Genera heatmaps (mapas de calor) del desempeño temporal del pronóstico
WRF-Chem contra observaciones SINAICA/INECC.

Ejes del heatmap
----------------
  Eje X (columnas) : meses del año (Ene – Dic) o semanas del año
  Eje Y (filas)    : ciudades del dominio

Cada celda muestra el valor de una métrica de evaluación, coloreada
por un gradiente que refleja el desempeño:

  Métricas disponibles (--metrica):
    POD   — Probabilidad de detección       [0..1]  verde=bueno
    FAR   — Tasa de falsas alarmas          [0..1]  verde=bajo (bueno)
    CSI   — Índice de éxito crítico         [0..1]  verde=bueno
    TSS   — Pierce Skill Score             [-1..1]  verde=bueno
    PC    — Porcentaje correcto             [0..1]  verde=bueno
    BIAS  — Sesgo de frecuencia            [0..∞]  verde=cerca de 1
    RMSE  — Raíz del error cuadrático medio  [0..∞]  verde=bajo
    R     — Correlación de Pearson         [-1..1]  verde=alto
    BIAS_cont — Sesgo continuo (mod − obs)   [−∞..∞]  verde=cerca de 0

Fuentes de datos
----------------
  Métricas dicotómicas (POD, FAR, CSI, TSS, PC, BIAS, BIAS_cont):
    Calculadas directamente desde combinado/ajustados/eval_*.csv

  Métricas continuas (RMSE, R):
    Calculadas directamente desde combinado/ajustados/eval_*.csv
    (no requieren el CSV de Taylor ni el CSV de auditoría dicotómica)

Temporadas climáticas marcadas en el eje X
------------------------------------------
  Secas frías     : nov – feb   (fondo gris azulado)
  Secas calientes : mar – may   (fondo amarillo claro)   ← máx O₃
  Lluvias         : jun – oct   (fondo verde claro)

Umbrales normativos utilizados (NOM-172-SEMARNAT-2023)
------------------------------------------------------
  Categoría "mala"     : O3=135ppb  PM10=132µg/m³  PM25=79µg/m³  SO2=185ppb
  Categoría "muy_mala" : O3=175ppb  PM10=213µg/m³  PM25=130µg/m³ SO2=304ppb

Uso
---
  # Heatmap mensual POD, todos los contaminantes, todas las ciudades
  python3 heatmap_desempeno.py

  # Solo O3, horizonte +24h, categoría muy mala
  python3 heatmap_desempeno.py --cont O3 --horizonte 24h --categoria muy_mala

  # Métrica RMSE (continua), resolución semanal
  python3 heatmap_desempeno.py --metrica RMSE --resolucion semana

  # Varias ciudades filtradas
  python3 heatmap_desempeno.py --ciudades CDMX Tula Pachuca

  # Rango de fechas
  python3 heatmap_desempeno.py --inicio 2026-01 --fin 2026-12

  # Guardar en un directorio específico
  python3 heatmap_desempeno.py --salida resultados_heatmap/

  python3 heatmap_desempeno.py --help

Salidas
-------
  heatmap_<CONT>_<METRICA>_<HOR>_<RESOLUCION>_<CATEGORIA>.png
  Una imagen PNG por cada contaminante procesado.

Dependencias
------------
  pip install pandas numpy matplotlib scipy

Autor  : Pipeline ddsinaica / WRF-Chem — ICAyCC, UNAM
Versión: 1.0.0 (2026-07)
"""

import argparse
import glob
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────

# Umbrales NOM-172-SEMARNAT-2023
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
    "O3":   {"nombre": "Ozono (O₃)",                      "unidad": "ppbv"},
    "PM10": {"nombre": "PM10",                             "unidad": "µg/m³"},
    "PM25": {"nombre": "PM2.5",                            "unidad": "µg/m³"},
    "SO2":  {"nombre": "Dióxido de azufre (SO₂)",         "unidad": "ppbv"},
}

HORIZONTES_MAP = {
    "24h": "mod_dia1", "+24h": "mod_dia1", "dia1": "mod_dia1",
    "48h": "mod_dia2", "+48h": "mod_dia2", "dia2": "mod_dia2",
    "72h": "mod_dia3", "+72h": "mod_dia3", "dia3": "mod_dia3",
}
HORIZONTES_LABEL = {
    "mod_dia1": "+24 h", "mod_dia2": "+48 h", "mod_dia3": "+72 h",
}

# Control de calidad (igual que informe_dicotomico.py y taylor_mensual.py)
LIMITES_VALIDOS = {
    "O3":   (0.0,  300.0,  0.0,  300.0),
    "PM10": (0.0, 1000.0,  0.0, 1000.0),
    "PM25": (0.0,  500.0,  0.0,  500.0),
    "SO2":  (0.0, 1000.0,  0.0, 1000.0),
}
CONTAMINANTES_FILTRO_IQR = {"PM25"}
IQR_FACTOR    = 3.0
UMBRAL_PM25   = 500.0   # techo absoluto obs PM2.5; puede sobreescribirse con --umbral-pm25
MIN_DIAS      = 5       # mínimo de pares válidos para calcular una celda

# Temporadas climáticas para el centro de México
TEMPORADAS = {
    "Secas frías":      {"meses": {11, 12, 1, 2},    "color": "#DDEEFF", "alpha": 0.35},
    "Secas calientes":  {"meses": {3, 4, 5},          "color": "#FFF9C4", "alpha": 0.40},
    "Lluvias":          {"meses": {6, 7, 8, 9, 10},  "color": "#DCEDC8", "alpha": 0.35},
}

DPI         = 150
NOMBRE_MES  = {
    1: "Ene", 2: "Feb", 3: "Mar", 4: "Abr", 5: "May", 6: "Jun",
    7: "Jul", 8: "Ago", 9: "Sep", 10: "Oct", 11: "Nov", 12: "Dic",
}

# ── Configuración por métrica ─────────────────────────────────────────────────
# cmap     : mapa de colores matplotlib
# vmin/vmax: rango de la escala de colores
# ideal    : valor perfecto (anotado en la barra de color)
# invert   : si True, el rojo es el valor alto (FAR, RMSE)
# fmt      : formato de las anotaciones numéricas en las celdas
METRICAS_CFG = {
    "POD":       {"cmap": "RdYlGn",        "vmin": 0.0,  "vmax": 1.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "FAR":       {"cmap": "RdYlGn_r",      "vmin": 0.0,  "vmax": 1.0,  "ideal": 0.0,  "invert": True,  "fmt": ".2f"},
    "CSI":       {"cmap": "RdYlGn",        "vmin": 0.0,  "vmax": 1.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "TSS":       {"cmap": "RdYlGn",        "vmin": -1.0, "vmax": 1.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "PC":        {"cmap": "RdYlGn",        "vmin": 0.0,  "vmax": 1.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "BIAS":      {"cmap": "RdYlGn",        "vmin": 0.0,  "vmax": 2.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "RMSE":      {"cmap": "RdYlGn_r",      "vmin": None, "vmax": None, "ideal": 0.0,  "invert": True,  "fmt": ".1f"},
    "R":         {"cmap": "RdYlGn",        "vmin": -1.0, "vmax": 1.0,  "ideal": 1.0,  "invert": False, "fmt": ".2f"},
    "BIAS_cont": {"cmap": "RdBu_r",        "vmin": None, "vmax": None, "ideal": 0.0,  "invert": False, "fmt": ".1f"},
}


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


# ──────────────────────────────────────────────────────────────────────────────
# CONTROL DE CALIDAD (igual que informe_dicotomico.py / taylor_mensual.py)
# ──────────────────────────────────────────────────────────────────────────────

def limpiar_serie(
    obs: np.ndarray,
    mod: np.ndarray,
    contaminante: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Elimina NaN, valores fuera de límites físicos y outliers IQR en PM2.5."""
    # 1. NaN / Inf
    mask = np.isfinite(obs) & np.isfinite(mod)
    obs, mod = obs[mask], mod[mask]

    # 2. Límites físicos
    lim = LIMITES_VALIDOS.get(contaminante)
    if lim:
        min_obs, max_obs, min_mod, max_mod = lim
        # Actualizar techo PM2.5 según CLI
        if contaminante == "PM25":
            max_obs = UMBRAL_PM25
        mask_lim = (
            (obs >= min_obs) & (obs <= max_obs) &
            (mod >= min_mod) & (mod <= max_mod)
        )
        obs, mod = obs[mask_lim], mod[mask_lim]

    # 3. Filtro IQR (solo PM2.5)
    if contaminante in CONTAMINANTES_FILTRO_IQR and len(obs) >= 4:
        q1, q3 = np.percentile(obs, [25, 75])
        iqr = q3 - q1
        if iqr > 0:
            mask_iqr = (obs >= q1 - IQR_FACTOR * iqr) & (obs <= q3 + IQR_FACTOR * iqr)
            obs, mod = obs[mask_iqr], mod[mask_iqr]

    return obs, mod


# ──────────────────────────────────────────────────────────────────────────────
# CÁLCULO DE MÉTRICAS
# ──────────────────────────────────────────────────────────────────────────────

def calcular_metrica(
    obs: np.ndarray,
    mod: np.ndarray,
    metrica: str,
    umbral: float,
) -> Optional[float]:
    """
    Calcula una única métrica escalar a partir de dos arrays obs/mod ya limpios.
    Retorna None si no hay suficientes datos.
    """
    if len(obs) < MIN_DIAS:
        return None

    def safe(num, den): return num / den if den > 0 else np.nan

    if metrica in ("POD", "FAR", "CSI", "TSS", "PC", "BIAS"):
        obs_ev = obs >= umbral
        mod_ev = mod >= umbral
        H = int(np.sum( obs_ev &  mod_ev))
        M = int(np.sum( obs_ev & ~mod_ev))
        F = int(np.sum(~obs_ev &  mod_ev))
        C = int(np.sum(~obs_ev & ~mod_ev))
        N = H + M + F + C

        if metrica == "POD":  v = safe(H, H + M)
        elif metrica == "FAR": v = safe(F, H + F)
        elif metrica == "CSI": v = safe(H, H + M + F)
        elif metrica == "TSS": v = safe(H, H + M) - safe(F, F + C)
        elif metrica == "PC":  v = safe(H + C, N)
        elif metrica == "BIAS": v = safe(H + F, H + M)

    elif metrica == "RMSE":
        v = float(np.sqrt(np.mean((mod - obs) ** 2)))

    elif metrica == "R":
        if np.std(obs) < 1e-10 or np.std(mod) < 1e-10:
            return None
        v, _ = stats.pearsonr(obs, mod)

    elif metrica == "BIAS_cont":
        v = float(np.mean(mod - obs))

    else:
        raise ValueError(f"Métrica no reconocida: {metrica}")

    return float(v) if np.isfinite(v) else None


# ──────────────────────────────────────────────────────────────────────────────
# LECTURA Y AGRUPACIÓN DE DATOS
# ──────────────────────────────────────────────────────────────────────────────

def leer_datos(
    directorio: str,
    ciudades_filtro: Optional[List[str]],
    inicio: Optional[str],
    fin: Optional[str],
    contaminante_filtro: Optional[str],
) -> pd.DataFrame:
    """
    Lee todos los eval_*.csv del directorio y devuelve un DataFrame largo con
    columnas: fecha, ciudad, contaminante, max_obs, mod_dia1, mod_dia2, mod_dia3.
    Aplica los filtros de ciudad, rango de fechas y contaminante.
    """
    archivos = sorted(glob.glob(os.path.join(directorio, "eval_*.csv")))
    if not archivos:
        sys.exit(f"[ERROR] No se encontraron eval_*.csv en '{directorio}'.")

    partes = []
    n_ok   = 0

    for ruta in archivos:
        meta = parsear_nombre(ruta)
        if meta is None:
            continue
        cont, ciudad = meta

        if cont not in CATEGORIAS_NOM172["mala"]:       # solo contaminantes conocidos
            continue
        if contaminante_filtro and cont != contaminante_filtro:
            continue
        if ciudades_filtro and ciudad not in ciudades_filtro:
            continue

        try:
            df = pd.read_csv(ruta, parse_dates=["Fecha"])
        except Exception:
            continue

        cols_req = {"Fecha", "max_obs", "mod_dia1", "mod_dia2", "mod_dia3"}
        if not cols_req.issubset(df.columns):
            continue

        df["ciudad"]       = ciudad
        df["contaminante"] = cont
        partes.append(df[["Fecha", "ciudad", "contaminante",
                           "max_obs", "mod_dia1", "mod_dia2", "mod_dia3"]])
        n_ok += 1

    if not partes:
        sys.exit("[ERROR] Ningún archivo coincidió con los filtros.")

    datos = pd.concat(partes, ignore_index=True)
    datos = datos.rename(columns={"Fecha": "fecha"})
    datos["fecha"] = pd.to_datetime(datos["fecha"])

    # Filtro de rango de fechas
    if inicio:
        datos = datos[datos["fecha"] >= pd.Timestamp(inicio + "-01")]
    if fin:
        datos = datos[datos["fecha"] <= pd.Timestamp(fin + "-01") + pd.offsets.MonthEnd(0)]

    print(f"[INFO] Registros cargados: {len(datos):,}  ({n_ok} archivos, "
          f"{datos['fecha'].dt.date.nunique()} días únicos)")
    return datos


# ──────────────────────────────────────────────────────────────────────────────
# CONSTRUCCIÓN DE LA MATRIZ DEL HEATMAP
# ──────────────────────────────────────────────────────────────────────────────

def construir_matriz_mensual(
    datos: pd.DataFrame,
    contaminante: str,
    col_mod: str,
    metrica: str,
    umbral: float,
    ciudades: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Construye un DataFrame (ciudades × YYYY-MM) con el valor de la métrica
    calculada mes a mes.  Retorna también la lista ordenada de etiquetas de mes.
    """
    datos_cont = datos[datos["contaminante"] == contaminante].copy()
    datos_cont["mes"] = datos_cont["fecha"].dt.to_period("M").astype(str)

    meses_disponibles = sorted(datos_cont["mes"].unique())
    matriz = {}

    for ciudad in ciudades:
        fila = {}
        df_c = datos_cont[datos_cont["ciudad"] == ciudad]
        for mes in meses_disponibles:
            df_m = df_c[df_c["mes"] == mes]
            if df_m.empty:
                fila[mes] = np.nan
                continue
            obs = np.asarray(df_m["max_obs"].values, float)
            mod = np.asarray(df_m[col_mod].values, float)
            obs, mod = limpiar_serie(obs, mod, contaminante)
            fila[mes] = calcular_metrica(obs, mod, metrica, umbral)
        matriz[ciudad] = fila

    df_mat = pd.DataFrame(matriz).T           # índice=ciudades, columnas=meses
    df_mat = df_mat.reindex(ciudades)         # orden canónico
    df_mat = df_mat[meses_disponibles]
    return df_mat, meses_disponibles


def construir_matriz_semanal(
    datos: pd.DataFrame,
    contaminante: str,
    col_mod: str,
    metrica: str,
    umbral: float,
    ciudades: List[str],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Construye la matriz (ciudades × semana-ISO) agrupando por semana ISO.
    """
    datos_cont = datos[datos["contaminante"] == contaminante].copy()
    datos_cont["semana"] = datos_cont["fecha"].dt.strftime("%G-W%V")   # ISO year-week

    semanas = sorted(datos_cont["semana"].unique())
    matriz  = {}

    for ciudad in ciudades:
        fila = {}
        df_c = datos_cont[datos_cont["ciudad"] == ciudad]
        for sem in semanas:
            df_s = df_c[df_c["semana"] == sem]
            if df_s.empty:
                fila[sem] = np.nan
                continue
            obs = np.asarray(df_s["max_obs"].values, float)
            mod = np.asarray(df_s[col_mod].values, float)
            obs, mod = limpiar_serie(obs, mod, contaminante)
            fila[sem] = calcular_metrica(obs, mod, metrica, umbral)
        matriz[ciudad] = fila

    df_mat = pd.DataFrame(matriz).T
    df_mat = df_mat.reindex(ciudades)
    df_mat = df_mat[semanas]
    return df_mat, semanas


# ──────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN
# ──────────────────────────────────────────────────────────────────────────────

def _etiquetas_meses(meses: List[str]) -> List[str]:
    """Convierte "2026-03" → "Mar 26"."""
    resultado = []
    for m in meses:
        anio, mm = m.split("-")
        resultado.append(f"{NOMBRE_MES[int(mm)]}\n{anio[2:]}")
    return resultado

def _etiquetas_semanas(semanas: List[str]) -> List[str]:
    """Muestra solo semanas múltiplo de 4 para no saturar el eje."""
    resultado = []
    for i, s in enumerate(semanas):
        resultado.append(s if i % 4 == 0 else "")
    return resultado

def _fondo_temporadas_mensual(ax: plt.Axes, meses: List[str]):
    """Pinta bandas de color de fondo según la temporada climática."""
    for i, mes in enumerate(meses):
        num_mes = int(mes.split("-")[1])
        for nombre_temp, cfg in TEMPORADAS.items():
            if num_mes in cfg["meses"]:
                ax.axvspan(i - 0.5, i + 0.5,
                           color=cfg["color"], alpha=cfg["alpha"], zorder=0)
                break

def _fondo_temporadas_semanal(ax: plt.Axes, semanas: List[str]):
    """Bandas de temporada para resolución semanal."""
    for i, sem in enumerate(semanas):
        try:
            # Calcular número de mes desde la semana ISO
            fecha = pd.Timestamp.fromisocalendar(
                int(sem[:4]), int(sem[6:]), 1
            )
            num_mes = fecha.month
        except Exception:
            continue
        for cfg in TEMPORADAS.values():
            if num_mes in cfg["meses"]:
                ax.axvspan(i - 0.5, i + 0.5,
                           color=cfg["color"], alpha=cfg["alpha"], zorder=0)
                break

def _leyenda_temporadas(ax: plt.Axes):
    """Añade parches de leyenda para las temporadas."""
    parches = [
        mpatches.Patch(color=cfg["color"], alpha=0.7, label=nombre)
        for nombre, cfg in TEMPORADAS.items()
    ]
    ax.legend(handles=parches, loc="upper left", bbox_to_anchor=(0, -0.14),
              ncol=3, fontsize=7.5, frameon=True, edgecolor="gray")


def _valor_a_color_bias(val: float, cfg: dict) -> str:
    """Para BIAS: mapea distancia a 1 en color de texto."""
    if val is None or np.isnan(val):
        return "gray"
    return "black"


def graficar_heatmap(
    df_mat: pd.DataFrame,
    etiquetas_x: List[str],
    columnas_x: List[str],
    contaminante: str,
    metrica: str,
    horizonte_label: str,
    resolucion: str,
    cat_key: str,
    umbral: float,
    dir_salida: str,
    anotar: bool = True,
):
    """
    Dibuja y guarda el heatmap.

    df_mat     : DataFrame (ciudades × etiquetas_x) con valores de la métrica
    etiquetas_x: etiquetas para el eje X
    columnas_x : claves internas (para ordenar)
    """
    cfg = METRICAS_CFG[metrica]
    meta = META_CONT[contaminante]
    n_cities = len(df_mat)
    n_cols   = len(columnas_x)

    # ── Tamaño de figura adaptativo ──────────────────────────────────────────
    ancho = max(10, n_cols * 0.85)
    alto  = max(3.5, n_cities * 0.72 + 2.5)
    fig, ax = plt.subplots(figsize=(ancho, alto))

    # ── Datos como array 2D ───────────────────────────────────────────────────
    data_arr = df_mat[columnas_x].values.astype(float)

    # Para BIAS: centrar la escala en 1
    vmin = cfg["vmin"]
    vmax = cfg["vmax"]
    if metrica == "BIAS":
        rango = np.nanmax(np.abs(data_arr - 1.0)) * 1.1 or 1.0
        vmin, vmax = max(0, 1 - rango), 1 + rango
        cmap_bias = plt.cm.get_cmap("RdYlGn")
        # Normalización centrada en 1
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vcenter=1.0, vmin=vmin, vmax=vmax)
    elif metrica == "BIAS_cont":
        lim = np.nanpercentile(np.abs(data_arr), 95) * 1.1 or 1.0
        vmin, vmax = -lim, lim
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vcenter=0.0, vmin=vmin, vmax=vmax)
    elif vmin is None or vmax is None:
        p5  = np.nanpercentile(data_arr,  5)
        p95 = np.nanpercentile(data_arr, 95)
        vmin = p5  - abs(p5)  * 0.1
        vmax = p95 + abs(p95) * 0.1
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = plt.Normalize(vmin=vmin, vmax=vmax)

    # ── Fondo de temporadas ───────────────────────────────────────────────────
    if resolucion == "mes":
        _fondo_temporadas_mensual(ax, columnas_x)
    else:
        _fondo_temporadas_semanal(ax, columnas_x)

    # ── Heatmap principal ─────────────────────────────────────────────────────
    im = ax.imshow(
        data_arr,
        cmap=cfg["cmap"],
        norm=norm,
        aspect="auto",
        interpolation="nearest",
        zorder=1,
    )

    # ── Anotaciones numéricas en cada celda ───────────────────────────────────
    if anotar:
        fmt = cfg["fmt"]
        for i in range(n_cities):
            for j in range(n_cols):
                val = data_arr[i, j]
                if np.isnan(val):
                    texto = "N/D"
                    color_txt = "gray"
                else:
                    texto = f"{val:{fmt}}"
                    # Color del texto: blanco sobre fondos oscuros, negro sobre claros
                    rgba = im.cmap(norm(val))
                    luminancia = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
                    color_txt  = "white" if luminancia < 0.50 else "black"
                ax.text(j, i, texto, ha="center", va="center",
                        fontsize=7.0, color=color_txt, zorder=2,
                        fontweight="bold" if not np.isnan(val) else "normal")

    # ── Ejes ─────────────────────────────────────────────────────────────────
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(etiquetas_x, fontsize=7.5, rotation=0)
    ax.set_yticks(range(n_cities))
    ax.set_yticklabels(df_mat.index.tolist(), fontsize=9)

    # Líneas de rejilla entre celdas
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_cities, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8, zorder=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # ── Barra de color ────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    unidad = meta["unidad"] if metrica in ("RMSE", "BIAS_cont") else ""
    cbar.set_label(f"{metrica}  {unidad}".strip(), fontsize=8.5)
    if cfg["ideal"] is not None:
        cbar.ax.axhline(y=norm(cfg["ideal"]), color="black",
                        linewidth=1.5, linestyle="--")
        cbar.ax.text(1.6, norm(cfg["ideal"]), f"ideal={cfg['ideal']}",
                     va="center", fontsize=7, transform=cbar.ax.transData)

    # ── Títulos ───────────────────────────────────────────────────────────────
    cat_txt = "Mala (naranja)" if cat_key == "mala" else "Muy Mala (rojo)"
    ax.set_title(
        f"Desempeño WRF-Chem — {meta['nombre']}    "
        f"Métrica: {metrica}    Horizonte: {horizonte_label}    "
        f"Umbral NOM-172: {umbral} {meta['unidad']} ({cat_txt})",
        fontsize=10, fontweight="bold", pad=10,
    )
    resolucion_txt = "mensual" if resolucion == "mes" else "semanal"
    ax.set_xlabel(f"Período ({resolucion_txt})", fontsize=9, labelpad=28)
    ax.set_ylabel("Ciudad", fontsize=9)

    # ── Leyenda de temporadas ─────────────────────────────────────────────────
    _leyenda_temporadas(ax)

    # ── Nota al pie ───────────────────────────────────────────────────────────
    fig.text(
        0.01, 0.01,
        f"Fuente: WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | "
        f"ddsinaica / ICAyCC, UNAM | N/D = < {MIN_DIAS} días válidos",
        fontsize=6.5, style="italic", color="gray",
    )

    plt.tight_layout(rect=[0, 0.06, 1, 1])

    # ── Guardar ───────────────────────────────────────────────────────────────
    os.makedirs(dir_salida, exist_ok=True)
    nombre = (f"heatmap_{contaminante}_{metrica}_{horizonte_label.replace(' ', '')}"
              f"_{resolucion}_{cat_key}.png")
    ruta   = os.path.join(dir_salida, nombre)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# ARGUMENTOS CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="heatmap_desempeno.py",
        description=(
            "Genera heatmaps (mes × ciudad o semana × ciudad) del desempeño "
            "temporal del pronóstico WRF-Chem vs SINAICA/INECC. "
            "Resalta automáticamente las temporadas climáticas del centro de México."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--entrada", "-i", default="combinado/ajustados",
                   help="Directorio con eval_*.csv (default: combinado/ajustados)")
    p.add_argument("--salida",  "-o", default="resultados_heatmap",
                   help="Directorio de salida (default: resultados_heatmap)")
    p.add_argument("--cont", "-C",    default=None,
                   choices=CONTAMINANTES_ORDEN,
                   help="Contaminante a procesar. Por defecto procesa los cuatro.")
    p.add_argument("--metrica", "-M", default="POD",
                   choices=list(METRICAS_CFG.keys()),
                   help="Métrica a visualizar (default: POD)")
    p.add_argument("--horizonte", "-H", default="24h",
                   help="Horizonte de pronóstico: 24h, 48h o 72h (default: 24h)")
    p.add_argument("--resolucion", "-r", default="mes",
                   choices=["mes", "semana"],
                   help="Resolución temporal del eje X (default: mes)")
    p.add_argument("--categoria", "-k", default="mala",
                   help="Categoría NOM-172: 'mala' o 'muy_mala' (default: mala)")
    p.add_argument("--ciudades", "-c", nargs="+", default=None,
                   help=f"Ciudades a incluir. Catálogo: {', '.join(CIUDADES_DOMINIO)}")
    p.add_argument("--inicio", default=None, metavar="YYYY-MM",
                   help="Mes de inicio del período (ej: 2026-01)")
    p.add_argument("--fin",    default=None, metavar="YYYY-MM",
                   help="Mes de fin del período (ej: 2026-12)")
    p.add_argument("--umbral-pm25", type=float, default=UMBRAL_PM25, metavar="UG_M3",
                   help=f"Techo absoluto PM2.5 en µg/m³ (default: {UMBRAL_PM25})")
    p.add_argument("--iqr-factor",  type=float, default=IQR_FACTOR,  metavar="K",
                   help=f"Factor IQR para outliers PM2.5 (default: {IQR_FACTOR})")
    p.add_argument("--min-dias",    type=int,   default=MIN_DIAS,    metavar="N",
                   help=f"Mínimo de días válidos por celda (default: {MIN_DIAS})")
    p.add_argument("--sin-anotar",  action="store_true",
                   help="Omitir las anotaciones numéricas en las celdas")
    p.add_argument("--dpi",         type=int,   default=DPI,
                   help=f"Resolución de los PNG (default: {DPI})")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# PUNTO DE ENTRADA
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    # Resolver parámetros globales
    global UMBRAL_PM25, IQR_FACTOR, MIN_DIAS, DPI
    UMBRAL_PM25 = args.umbral_pm25
    IQR_FACTOR  = args.iqr_factor
    MIN_DIAS    = args.min_dias
    DPI         = args.dpi

    # Categoría
    cat_key = ALIAS_CATEGORIA.get(args.categoria.lower().replace(" ", "_"))
    if cat_key is None:
        sys.exit(f"[ERROR] Categoría no reconocida: '{args.categoria}'\n"
                 f"        Opciones: mala, muy_mala")

    # Horizonte
    col_mod = HORIZONTES_MAP.get(args.horizonte.lower().replace("+", ""))
    if col_mod is None:
        sys.exit(f"[ERROR] Horizonte no reconocido: '{args.horizonte}'\n"
                 f"        Opciones: 24h, 48h, 72h")
    horizonte_label = HORIZONTES_LABEL[col_mod]

    # Ciudades
    ciudades_filtro = parsear_lista_ciudades(args.ciudades)
    ciudades_activas = ciudades_filtro or CIUDADES_DOMINIO

    # Contaminantes a procesar
    conts = [args.cont] if args.cont else CONTAMINANTES_ORDEN

    print("=" * 66)
    print("  Heatmap de desempeño — WRF-Chem / SINAICA")
    print("=" * 66)
    print(f"  Entrada      : {args.entrada}")
    print(f"  Salida       : {args.salida}")
    print(f"  Métrica      : {args.metrica}")
    print(f"  Horizonte    : {horizonte_label}")
    print(f"  Resolución   : {args.resolucion}")
    print(f"  Categoría    : {cat_key} (NOM-172-SEMARNAT-2023)")
    print(f"  Ciudades     : {ciudades_activas}")
    print(f"  Período      : {args.inicio or 'inicio'} – {args.fin or 'fin'}")
    print(f"  QC PM2.5     : techo={UMBRAL_PM25} µg/m³  IQR·k={IQR_FACTOR}")
    print("=" * 66)

    # Leer datos
    datos = leer_datos(
        args.entrada, ciudades_filtro,
        args.inicio, args.fin,
        args.cont,
    )

    # Generar un heatmap por contaminante
    for cont in conts:
        if cont not in datos["contaminante"].unique():
            print(f"[WARN] {cont}: sin datos — omitido.")
            continue

        umbral = CATEGORIAS_NOM172[cat_key][cont]

        # Filtrar ciudades con datos para este contaminante
        ciudades_con_datos = [
            c for c in ciudades_activas
            if c in datos.loc[datos["contaminante"] == cont, "ciudad"].unique()
        ]
        if not ciudades_con_datos:
            print(f"[WARN] {cont}: ninguna ciudad con datos — omitido.")
            continue

        print(f"\n── {cont} (umbral {cat_key}: {umbral} {META_CONT[cont]['unidad']}) ──")

        if args.resolucion == "mes":
            df_mat, cols = construir_matriz_mensual(
                datos, cont, col_mod, args.metrica, umbral, ciudades_con_datos
            )
            etiquetas_x = _etiquetas_meses(cols)
        else:
            df_mat, cols = construir_matriz_semanal(
                datos, cont, col_mod, args.metrica, umbral, ciudades_con_datos
            )
            etiquetas_x = _etiquetas_semanas(cols)

        # Eliminar ciudades sin ningún dato
        df_mat = df_mat.dropna(how="all")
        if df_mat.empty:
            print(f"[WARN] {cont}: matriz vacía tras QC — omitido.")
            continue

        graficar_heatmap(
            df_mat       = df_mat,
            etiquetas_x  = etiquetas_x,
            columnas_x   = cols,
            contaminante = cont,
            metrica      = args.metrica,
            horizonte_label = horizonte_label,
            resolucion   = args.resolucion,
            cat_key      = cat_key,
            umbral       = umbral,
            dir_salida   = args.salida,
            anotar       = not args.sin_anotar,
        )

    print(f"\n[DONE] Heatmaps guardados en '{args.salida}/'")


if __name__ == "__main__":
    main()

