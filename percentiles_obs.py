#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
percentiles_obs.py
==================
Calcula y visualiza los percentiles empíricos de las observaciones SINAICA
y los compara contra los percentiles del modelo WRF-Chem por ciudad, mes y
contaminante.

Motivación científica
---------------------
Los umbrales normativos (NOM-172-SEMARNAT-2023) son nacionales y no reflejan
la climatología local de concentración. El P90 de O₃ en Tula puede estar 40 ppb
por encima del P90 de Cuernavaca. Usar el mismo umbral para evaluar el modelo en
ambas ciudades mezcla distribuciones muy distintas.

Este módulo responde tres preguntas:
    1. ¿Cuál es la distribución empírica de concentraciones por ciudad y mes?
       → permite definir umbrales locales adaptativos (ver curva_habilidad.py)
    2. ¿El modelo reproduce los percentiles altos de la distribución observada?
       → un R alto puede coexistir con subestimación sistemática del P90
    3. ¿Hay diferencias estacionales (secas vs lluvias) en la distribución?

Salidas
-------
    percentiles_<CONT>_heatmap_obs.png
        Heatmap (ciudad × mes) del percentil seleccionado de las observaciones.
        Revela la climatología local de concentración y la variación estacional.

    percentiles_<CONT>_heatmap_bias.png
        Heatmap (ciudad × mes) del sesgo percentílico: P_mod / P_obs − 1.
        Verde = modelo reproduce bien; rojo = subestima; azul = sobreestima.

    percentiles_<CONT>_scatter_QQ.png
        Diagrama Q-Q (cuantil obs vs cuantil mod) con una curva por temporada.
        La línea 1:1 representa reproducción perfecta de la distribución.

    percentiles_<CONT>_serie_mensual.png
        Serie temporal de los percentiles P25/P50/P75/P90 para obs y mod,
        con bandas sombreadas entre P25-P75 (rango intercuartílico).

    tabla_percentiles_<CONT>.csv
        ciudad, mes, temporada, contaminante, horizonte,
        N, P10_obs, P25_obs, P50_obs, P75_obs, P90_obs, P95_obs, P99_obs,
        P10_mod, P25_mod, P50_mod, P75_mod, P90_mod, P95_mod, P99_mod,
        ratio_P50, ratio_P90, ratio_IQR, sesgo_P90 (P_mod/P_obs − 1)

Uso
---
    # Todos los contaminantes, todos los horizontes
    python3 percentiles_obs.py

    # Solo O3, horizonte +24h
    python3 percentiles_obs.py --cont O3 --horizonte 24h

    # Solo ciertas ciudades
    python3 percentiles_obs.py --ciudades CDMX Tula Pachuca

    # Percentil de referencia para el heatmap de sesgo (default: P90)
    python3 percentiles_obs.py --percentil 95

    # Período específico
    python3 percentiles_obs.py --inicio 2026-01 --fin 2026-12

    python3 percentiles_obs.py --help

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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────

# Umbrales NOM-172-SEMARNAT-2023 — se marcan en las gráficas como referencia
NOM172 = {
    "mala":     {"O3": 135.0, "PM10": 132.0, "PM25":  79.0, "SO2": 185.0},
    "muy_mala": {"O3": 175.0, "PM10": 213.0, "PM25": 130.0, "SO2": 304.0},
}

META_CONT = {
    "O3":   {"nombre": "Ozono (O₃)",              "unidad": "ppbv",  "color": "#1f77b4"},
    "PM10": {"nombre": "PM10",                    "unidad": "µg/m³", "color": "#d62728"},
    "PM25": {"nombre": "PM2.5",                   "unidad": "µg/m³", "color": "#ff7f0e"},
    "SO2":  {"nombre": "Dióxido de azufre (SO₂)", "unidad": "ppbv",  "color": "#9467bd"},
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

HORIZONTES_MAP = {
    "24h": "mod_dia1", "+24h": "mod_dia1", "dia1": "mod_dia1",
    "48h": "mod_dia2", "+48h": "mod_dia2", "dia2": "mod_dia2",
    "72h": "mod_dia3", "+72h": "mod_dia3", "dia3": "mod_dia3",
    "todos": None,
}
HORIZONTES_LABELS = {
    "mod_dia1": "+24 h", "mod_dia2": "+48 h", "mod_dia3": "+72 h",
}
COLORES_HOR = {
    "mod_dia1": "#1f77b4", "mod_dia2": "#ff7f0e", "mod_dia3": "#d62728",
}

# Temporadas climáticas del centro de México
TEMPORADAS = {
    "Secas frías":     {"meses": {11, 12, 1, 2}, "color": "#4A90D9", "marker": "o"},
    "Secas calientes": {"meses": {3, 4, 5},       "color": "#F5A623", "marker": "s"},
    "Lluvias":         {"meses": {6, 7, 8, 9, 10},"color": "#2d8a2d", "marker": "^"},
}

# Control de calidad (idéntico al resto del pipeline)
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

NOMBRE_MES = {
    1: "Ene", 2: "Feb",  3: "Mar",  4: "Abr",
    5: "May", 6: "Jun",  7: "Jul",  8: "Ago",
    9: "Sep", 10: "Oct", 11: "Nov", 12: "Dic",
}

# Percentiles calculados siempre
PCTS = [10, 25, 50, 75, 90, 95, 99]

DPI = 150


# ──────────────────────────────────────────────────────────────────────────────
# UTILIDADES
# ──────────────────────────────────────────────────────────────────────────────

def normalizar_ciudad(nombre: str) -> Optional[str]:
    return _CIUDADES_NORM.get(nombre.strip().lower())

def parsear_lista_ciudades(valor) -> Optional[List[str]]:
    if not valor: return None
    crudos = []
    for item in ([valor] if isinstance(valor, str) else valor):
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
    for nombre, cfg in TEMPORADAS.items():
        if mes in cfg["meses"]: return nombre
    return "Lluvias"

def etiqueta_mes(mes_str: str) -> str:
    """'2026-04' → 'Abr\n26'"""
    anio, mm = mes_str.split("-")
    return f"{NOMBRE_MES[int(mm)]}\n{anio[2:]}"


# ──────────────────────────────────────────────────────────────────────────────
# CONTROL DE CALIDAD
# ──────────────────────────────────────────────────────────────────────────────

def limpiar_serie(obs: np.ndarray, mod: np.ndarray,
                  cont: str) -> Tuple[np.ndarray, np.ndarray]:
    """NaN → límites físicos → filtro IQR sobre obs PM2.5."""
    mask = np.isfinite(obs) & np.isfinite(mod)
    obs, mod = obs[mask], mod[mask]
    lim = LIMITES_VALIDOS.get(cont)
    if lim:
        min_o, max_o, min_m, max_m = lim
        if cont == "PM25": max_o = UMBRAL_PM25
        mask2 = (obs>=min_o)&(obs<=max_o)&(mod>=min_m)&(mod<=max_m)
        obs, mod = obs[mask2], mod[mask2]
    if cont in CONTAMINANTES_FILTRO_IQR and len(obs) >= 4:
        q1, q3 = np.percentile(obs, [25, 75])
        iqr = q3 - q1
        if iqr > 0:
            mask3 = (obs >= q1-IQR_FACTOR*iqr) & (obs <= q3+IQR_FACTOR*iqr)
            obs, mod = obs[mask3], mod[mask3]
    return obs, mod


# ──────────────────────────────────────────────────────────────────────────────
# LECTURA DE DATOS
# ──────────────────────────────────────────────────────────────────────────────

def leer_datos(directorio: str, ciudades_filtro: Optional[List[str]],
               inicio: Optional[str], fin: Optional[str],
               cont_filtro: Optional[str]) -> pd.DataFrame:
    archivos = sorted(glob.glob(os.path.join(directorio, "eval_*.csv")))
    if not archivos:
        sys.exit(f"[ERROR] No se encontraron eval_*.csv en '{directorio}'.")
    partes = []
    for ruta in archivos:
        meta = parsear_nombre(ruta)
        if not meta: continue
        cont, ciudad = meta
        if cont not in META_CONT: continue
        if cont_filtro and cont != cont_filtro: continue
        if ciudades_filtro and ciudad not in ciudades_filtro: continue
        try:
            df = pd.read_csv(ruta, parse_dates=["Fecha"])
        except Exception:
            continue
        if not {"Fecha","max_obs","mod_dia1","mod_dia2","mod_dia3"}.issubset(df.columns):
            continue
        df["ciudad"] = ciudad; df["contaminante"] = cont
        partes.append(df[["Fecha","ciudad","contaminante",
                           "max_obs","mod_dia1","mod_dia2","mod_dia3"]])
    if not partes:
        sys.exit("[ERROR] Ningún archivo coincidió con los filtros.")
    datos = pd.concat(partes, ignore_index=True).rename(columns={"Fecha":"fecha"})
    datos["fecha"] = pd.to_datetime(datos["fecha"])
    if inicio:
        datos = datos[datos["fecha"] >= pd.Timestamp(inicio+"-01")]
    if fin:
        datos = datos[datos["fecha"] <= pd.Timestamp(fin+"-01")+pd.offsets.MonthEnd(0)]
    datos["mes"]       = datos["fecha"].dt.month
    datos["mes_str"]   = datos["fecha"].dt.to_period("M").astype(str)
    datos["temporada"] = datos["mes"].map(mes_a_temporada)
    print(f"[INFO] Registros: {len(datos):,}  "
          f"({datos['fecha'].dt.date.nunique()} días, "
          f"{datos['ciudad'].nunique()} ciudades)")
    return datos


# ──────────────────────────────────────────────────────────────────────────────
# CÁLCULO DE PERCENTILES
# ──────────────────────────────────────────────────────────────────────────────

def calcular_percentiles_serie(arr: np.ndarray) -> Dict:
    """Calcula los PCTS percentiles de un array. Retorna dict P10..P99."""
    if len(arr) < MIN_DIAS:
        return {f"P{p}": np.nan for p in PCTS}
    vals = np.percentile(arr, PCTS)
    return {f"P{p}": round(float(v), 3) for p, v in zip(PCTS, vals)}


def construir_tabla_percentiles(datos: pd.DataFrame,
                                 cont: str,
                                 cols_mod: List[str]) -> pd.DataFrame:
    """
    Construye el DataFrame de percentiles por ciudad × mes × horizonte.

    Columnas: ciudad, mes_str, mes, temporada, contaminante, horizonte, N,
              P10_obs..P99_obs, P10_mod..P99_mod,
              ratio_P50, ratio_P90, ratio_IQR, sesgo_P90
    """
    df_cont = datos[datos["contaminante"] == cont].copy()
    filas = []

    ciudades = sorted(df_cont["ciudad"].unique(),
                      key=lambda c: CIUDADES_DOMINIO.index(c)
                      if c in CIUDADES_DOMINIO else 99)
    meses    = sorted(df_cont["mes_str"].unique())

    for ciudad in ciudades:
        dfc = df_cont[df_cont["ciudad"] == ciudad]
        for mes_str in meses:
            dfm = dfc[dfc["mes_str"] == mes_str]
            if dfm.empty: continue
            mes_num = int(mes_str.split("-")[1])
            temporada = mes_a_temporada(mes_num)

            for col_mod in cols_mod:
                hor_label = HORIZONTES_LABELS[col_mod]
                obs = np.asarray(dfm["max_obs"].values, float)
                mod = np.asarray(dfm[col_mod].values, float)
                obs, mod = limpiar_serie(obs, mod, cont)

                if len(obs) < MIN_DIAS:
                    continue

                pct_obs = calcular_percentiles_serie(obs)
                pct_mod = calcular_percentiles_serie(mod)

                fila = {
                    "ciudad":       ciudad,
                    "mes_str":      mes_str,
                    "mes":          mes_num,
                    "temporada":    temporada,
                    "contaminante": cont,
                    "horizonte":    hor_label,
                    "N":            len(obs),
                }
                for p in PCTS:
                    fila[f"P{p}_obs"] = pct_obs[f"P{p}"]
                    fila[f"P{p}_mod"] = pct_mod[f"P{p}"]

                # Métricas derivadas
                p50_o = pct_obs["P50"]; p50_m = pct_mod["P50"]
                p90_o = pct_obs["P90"]; p90_m = pct_mod["P90"]
                p25_o = pct_obs["P25"]; p75_o = pct_obs["P75"]
                p25_m = pct_mod["P25"]; p75_m = pct_mod["P75"]

                fila["ratio_P50"]  = round(p50_m/p50_o, 4) if p50_o > 0 else np.nan
                fila["ratio_P90"]  = round(p90_m/p90_o, 4) if p90_o > 0 else np.nan
                fila["sesgo_P90"]  = round(p90_m/p90_o - 1, 4) if p90_o > 0 else np.nan
                iqr_obs = p75_o - p25_o
                iqr_mod = p75_m - p25_m
                fila["ratio_IQR"]  = round(iqr_mod/iqr_obs, 4) if iqr_obs > 0 else np.nan

                filas.append(fila)

    df_out = pd.DataFrame(filas)
    print(f"[INFO] {cont}: {len(df_out)} combinaciones ciudad×mes×horizonte")
    return df_out


# ──────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN — HEATMAP DE PERCENTIL OBSERVADO
# ──────────────────────────────────────────────────────────────────────────────

def graficar_heatmap_obs(df_tabla: pd.DataFrame, cont: str,
                          pct_ref: int, col_mod: str,
                          dir_salida: str):
    """
    Heatmap (ciudad × mes) del percentil P_pct_ref de las observaciones.
    Cada celda muestra el valor observado; las celdas del modelo se superponen
    como texto secundario más pequeño.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h      = df_tabla[df_tabla["horizonte"] == hor_label].copy()
    if df_h.empty:
        print(f"[WARN] {cont}/{hor_label}: sin datos para heatmap obs.")
        return

    ciudades  = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses     = sorted(df_h["mes_str"].unique())
    n_c, n_m  = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    # Matrices: filas=ciudades, cols=meses
    mat_obs  = np.full((n_c, n_m), np.nan)
    mat_mod  = np.full((n_c, n_m), np.nan)

    for i, ciudad in enumerate(ciudades):
        for j, mes_str in enumerate(meses):
            row = df_h[(df_h["ciudad"]==ciudad) & (df_h["mes_str"]==mes_str)]
            if row.empty: continue
            mat_obs[i, j] = row[f"P{pct_ref}_obs"].values[0]
            mat_mod[i, j] = row[f"P{pct_ref}_mod"].values[0]

    fig, ax = plt.subplots(figsize=(max(8, n_m*0.9), max(4, n_c*0.75 + 2)))

    # Fondo de temporadas en el eje X
    for j, mes_str in enumerate(meses):
        mm = int(mes_str.split("-")[1])
        for nombre_t, cfg_t in TEMPORADAS.items():
            if mm in cfg_t["meses"]:
                ax.axvspan(j-0.5, j+0.5, color=cfg_t["color"],
                           alpha=0.18, zorder=0)
                break

    vmin = np.nanmin(mat_obs); vmax = np.nanmax(mat_obs)
    im = ax.imshow(mat_obs, cmap="YlOrRd", aspect="auto",
                   vmin=vmin, vmax=vmax, interpolation="nearest", zorder=1)

    # Anotaciones: obs (grande) + mod (pequeña, gris)
    for i in range(n_c):
        for j in range(n_m):
            v_obs = mat_obs[i, j]
            v_mod = mat_mod[i, j]
            if np.isnan(v_obs):
                ax.text(j, i, "N/D", ha="center", va="center",
                        fontsize=7, color="gray", zorder=3)
                continue
            rgba = im.cmap(im.norm(v_obs))
            lum  = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            col_txt = "white" if lum < 0.50 else "black"
            ax.text(j, i-0.12, f"{v_obs:.1f}", ha="center", va="center",
                    fontsize=8.5, fontweight="bold", color=col_txt, zorder=3)
            if not np.isnan(v_mod):
                ax.text(j, i+0.25, f"mod:{v_mod:.1f}", ha="center", va="center",
                        fontsize=6.5, color=col_txt, alpha=0.80, zorder=3)

    # Rejilla
    ax.set_xticks(np.arange(-0.5, n_m, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0, zorder=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_xticks(range(n_m))
    ax.set_xticklabels([etiqueta_mes(m) for m in meses], fontsize=8)
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(ciudades, fontsize=9)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(f"P{pct_ref} observado  [{meta['unidad']}]", fontsize=9)

    # Marcar umbrales NOM-172 en la barra de color
    for cat_key, cat_vals in NOM172.items():
        u_nom = cat_vals.get(cont)
        if u_nom is None or u_nom < vmin or u_nom > vmax: continue
        col_nom = "#FF7E00" if cat_key=="mala" else "#CC0000"
        cbar.ax.axhline(y=(u_nom-vmin)/(vmax-vmin),
                        color=col_nom, lw=1.5, ls="--")
        cbar.ax.text(1.6, (u_nom-vmin)/(vmax-vmin),
                     f"NOM {cat_key.replace('_',' ')} ({u_nom})",
                     va="center", fontsize=6.5, color=col_nom,
                     transform=cbar.ax.transData)

    # Leyenda de temporadas
    parches = [mpatches.Patch(color=cfg["color"], alpha=0.6, label=nombre)
               for nombre, cfg in TEMPORADAS.items()]
    ax.legend(handles=parches, loc="upper left",
              bbox_to_anchor=(0, -0.12), ncol=3, fontsize=7.5, frameon=True)

    ax.set_title(
        f"P{pct_ref} observado y modelado — {meta['nombre']}\n"
        f"Horizonte {hor_label}  |  "
        f"(grande=obs, pequeño=mod:{hor_label})",
        fontsize=10, fontweight="bold", pad=10,
    )
    ax.set_xlabel("Período mensual", fontsize=9, labelpad=26)
    ax.set_ylabel("Ciudad", fontsize=9)
    fig.text(0.01, 0.01,
             f"Fuente: WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | "
             f"ICAyCC, UNAM | N/D = < {MIN_DIAS} días",
             fontsize=6.5, style="italic", color="gray")

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    nombre = f"percentiles_{cont}_P{pct_ref}_heatmap_obs_{hor_label.replace(' ','')}.png"
    ruta   = os.path.join(dir_salida, nombre)
    os.makedirs(dir_salida, exist_ok=True)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN — HEATMAP DE SESGO PERCENTÍLICO
# ──────────────────────────────────────────────────────────────────────────────

def graficar_heatmap_sesgo(df_tabla: pd.DataFrame, cont: str,
                            pct_ref: int, col_mod: str,
                            dir_salida: str):
    """
    Heatmap (ciudad × mes) del sesgo percentílico: P_mod/P_obs − 1.
    Divergente: rojo = subestimación (mod < obs); azul = sobreestimación.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h      = df_tabla[df_tabla["horizonte"] == hor_label].copy()
    if df_h.empty: return

    ciudades = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses    = sorted(df_h["mes_str"].unique())
    n_c, n_m = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    mat_sesgo = np.full((n_c, n_m), np.nan)
    mat_n     = np.full((n_c, n_m), 0, dtype=int)

    for i, ciudad in enumerate(ciudades):
        for j, mes_str in enumerate(meses):
            row = df_h[(df_h["ciudad"]==ciudad) & (df_h["mes_str"]==mes_str)]
            if row.empty: continue
            mat_sesgo[i, j] = row["sesgo_P90"].values[0] if pct_ref == 90 \
                              else (row[f"P{pct_ref}_mod"].values[0] /
                                    row[f"P{pct_ref}_obs"].values[0] - 1
                                    if row[f"P{pct_ref}_obs"].values[0] > 0
                                    else np.nan)
            mat_n[i, j] = row["N"].values[0]

    # Escala simétrica centrada en 0
    lim = np.nanpercentile(np.abs(mat_sesgo[np.isfinite(mat_sesgo)]), 95) if \
          np.any(np.isfinite(mat_sesgo)) else 0.5
    lim = max(lim, 0.05)

    fig, ax = plt.subplots(figsize=(max(8, n_m*0.9), max(4, n_c*0.75 + 2)))

    for j, mes_str in enumerate(meses):
        mm = int(mes_str.split("-")[1])
        for nombre_t, cfg_t in TEMPORADAS.items():
            if mm in cfg_t["meses"]:
                ax.axvspan(j-0.5, j+0.5, color=cfg_t["color"],
                           alpha=0.15, zorder=0)
                break

    from matplotlib.colors import TwoSlopeNorm
    norm_div = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
    im = ax.imshow(mat_sesgo, cmap="RdBu", aspect="auto",
                   norm=norm_div, interpolation="nearest", zorder=1)

    for i in range(n_c):
        for j in range(n_m):
            v = mat_sesgo[i, j]
            n = mat_n[i, j]
            if np.isnan(v):
                ax.text(j, i, "N/D", ha="center", va="center",
                        fontsize=7, color="gray", zorder=3)
                continue
            rgba = im.cmap(norm_div(v))
            lum  = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            col_txt = "white" if lum < 0.45 else "black"
            signo = "+" if v >= 0 else ""
            ax.text(j, i-0.12, f"{signo}{v*100:.1f}%",
                    ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color=col_txt, zorder=3)
            ax.text(j, i+0.25, f"n={n}", ha="center", va="center",
                    fontsize=6, color=col_txt, alpha=0.75, zorder=3)

    ax.set_xticks(np.arange(-0.5, n_m, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_c, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0, zorder=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xticks(range(n_m))
    ax.set_xticklabels([etiqueta_mes(m) for m in meses], fontsize=8)
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(ciudades, fontsize=9)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(f"Sesgo P{pct_ref}  (mod/obs − 1)", fontsize=9)
    cbar.ax.axhline(y=norm_div(0), color="black", lw=1.5, ls="--")

    parches = [mpatches.Patch(color=cfg["color"], alpha=0.6, label=nombre)
               for nombre, cfg in TEMPORADAS.items()]
    ax.legend(handles=parches, loc="upper left",
              bbox_to_anchor=(0, -0.12), ncol=3, fontsize=7.5)

    ax.set_title(
        f"Sesgo percentílico P{pct_ref} — {meta['nombre']}\n"
        f"Horizonte {hor_label}  |  "
        f"Rojo = modelo subestima  ·  Azul = modelo sobreestima",
        fontsize=10, fontweight="bold", pad=10,
    )
    ax.set_xlabel("Período mensual", fontsize=9, labelpad=26)
    ax.set_ylabel("Ciudad", fontsize=9)
    fig.text(0.01, 0.01,
             "WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | ICAyCC, UNAM",
             fontsize=6.5, style="italic", color="gray")

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    nombre = f"percentiles_{cont}_P{pct_ref}_heatmap_sesgo_{hor_label.replace(' ','')}.png"
    ruta   = os.path.join(dir_salida, nombre)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN — DIAGRAMA Q-Q POR TEMPORADA
# ──────────────────────────────────────────────────────────────────────────────

def graficar_qq_temporada(datos: pd.DataFrame, cont: str,
                           col_mod: str, dir_salida: str):
    """
    Diagrama Q-Q (cuantil obs vs cuantil mod) con una curva por temporada.
    Todos los días y ciudades se acumulan por temporada.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_cont   = datos[datos["contaminante"] == cont].copy()

    fig, ax = plt.subplots(figsize=(7, 6.5))

    # Línea 1:1
    all_vals = []
    for temp, cfg in TEMPORADAS.items():
        dft = df_cont[df_cont["temporada"] == temp]
        if dft.empty: continue
        obs = np.asarray(dft["max_obs"].values, float)
        mod = np.asarray(dft[col_mod].values, float)
        obs, mod = limpiar_serie(obs, mod, cont)
        if len(obs) < MIN_DIAS: continue
        all_vals.extend(obs.tolist() + mod.tolist())

    if not all_vals:
        plt.close(fig); return

    vmin = np.percentile(all_vals, 1)
    vmax = np.percentile(all_vals, 99)
    ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1.5, alpha=0.7,
            label="Línea 1:1 (reproducción perfecta)", zorder=2)

    # Banda ±15 % alrededor de la línea 1:1
    ax.fill_between([vmin, vmax],
                    [vmin*0.85, vmax*0.85],
                    [vmin*1.15, vmax*1.15],
                    color="gray", alpha=0.10, label="±15%", zorder=1)

    # Curva Q-Q por temporada
    for temp, cfg in TEMPORADAS.items():
        dft = df_cont[df_cont["temporada"] == temp]
        if dft.empty: continue
        obs = np.asarray(dft["max_obs"].values, float)
        mod = np.asarray(dft[col_mod].values, float)
        obs, mod = limpiar_serie(obs, mod, cont)
        if len(obs) < MIN_DIAS: continue

        # Percentiles para el diagrama Q-Q
        pcts_qq = np.linspace(1, 99, 99)
        q_obs   = np.percentile(obs, pcts_qq)
        q_mod   = np.percentile(mod, pcts_qq)

        ax.plot(q_obs, q_mod, color=cfg["color"], lw=2.0,
                label=f"{temp}  (n={len(obs)})", zorder=3)

        # Marcar percentiles clave
        for p_mark, marker in [(50, "o"), (90, "D"), (95, "^")]:
            idx = int(p_mark) - 1
            ax.scatter(q_obs[idx], q_mod[idx],
                       color=cfg["color"], s=60,
                       marker=marker, edgecolors="white",
                       linewidths=0.8, zorder=5)
            ax.annotate(f"P{p_mark}",
                        xy=(q_obs[idx], q_mod[idx]),
                        xytext=(q_obs[idx]+vmax*0.01, q_mod[idx]+vmax*0.01),
                        fontsize=6.5, color=cfg["color"], zorder=6)

    # Marcar umbrales NOM-172 con líneas verticales
    for cat_key, cat_vals in NOM172.items():
        u_nom = cat_vals.get(cont)
        if u_nom is None or u_nom < vmin or u_nom > vmax: continue
        col_nom = "#FF7E00" if cat_key=="mala" else "#CC0000"
        ax.axvline(u_nom, color=col_nom, lw=1.2, ls="--", alpha=0.75,
                   label=f"NOM-172 {cat_key.replace('_',' ')} ({u_nom} {meta['unidad']})")
        ax.axhline(u_nom, color=col_nom, lw=1.2, ls=":", alpha=0.50)

    ax.set_xlim(vmin, vmax); ax.set_ylim(vmin, vmax)
    ax.set_xlabel(f"Cuantiles observados  [{meta['unidad']}]", fontsize=10)
    ax.set_ylabel(f"Cuantiles modelo {hor_label}  [{meta['unidad']}]", fontsize=10)
    ax.set_title(
        f"Diagrama Q-Q por temporada — {meta['nombre']}\n"
        f"Horizonte {hor_label}  |  Puntos: ●=P50  ◆=P90  ▲=P95",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.90)
    ax.grid(True, color="#DDD", lw=0.5)
    ax.tick_params(labelsize=8.5)
    ax.set_aspect("equal")
    fig.text(0.01, 0.01,
             "WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | ICAyCC, UNAM",
             fontsize=6.5, style="italic", color="gray")

    plt.tight_layout()
    nombre = f"percentiles_{cont}_QQ_{hor_label.replace(' ','')}.png"
    ruta   = os.path.join(dir_salida, nombre)
    os.makedirs(dir_salida, exist_ok=True)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN — SERIE TEMPORAL DE PERCENTILES
# ──────────────────────────────────────────────────────────────────────────────

def graficar_serie_percentiles(df_tabla: pd.DataFrame, cont: str,
                                col_mod: str, ciudad: str,
                                dir_salida: str):
    """
    Serie temporal mensual de P25/P50/P75/P90 para una ciudad.
    Obs: línea continua. Mod: línea punteada del mismo color.
    Banda sombreada: rango intercuartílico (P25–P75).
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_c      = df_tabla[(df_tabla["horizonte"]==hor_label) &
                          (df_tabla["ciudad"]==ciudad)].sort_values("mes_str")
    if df_c.empty: return

    meses_str = df_c["mes_str"].tolist()
    x         = np.arange(len(meses_str))
    etqs      = [etiqueta_mes(m) for m in meses_str]

    fig, ax = plt.subplots(figsize=(max(8, len(x)*0.9), 5.5))

    # Fondo de temporadas
    for j, mes_str in enumerate(meses_str):
        mm = int(mes_str.split("-")[1])
        for nombre_t, cfg_t in TEMPORADAS.items():
            if mm in cfg_t["meses"]:
                ax.axvspan(j-0.5, j+0.5, color=cfg_t["color"],
                           alpha=0.15, zorder=0)
                break

    colores_pct = {
        "P25": "#74b9ff", "P50": "#0984e3",
        "P75": "#fdcb6e", "P90": "#d63031",
    }

    for pct_lbl, col in colores_pct.items():
        obs_vals = df_c[f"{pct_lbl}_obs"].values.astype(float)
        mod_vals = df_c[f"{pct_lbl}_mod"].values.astype(float)

        ax.plot(x, obs_vals, color=col, lw=2.0, marker="o", ms=5,
                label=f"{pct_lbl} obs", zorder=4)
        ax.plot(x, mod_vals, color=col, lw=1.5, ls="--", marker="s", ms=4,
                alpha=0.80, label=f"{pct_lbl} mod {hor_label}", zorder=3)

    # Banda IQR observado (P25-P75)
    p25_obs = df_c["P25_obs"].values.astype(float)
    p75_obs = df_c["P75_obs"].values.astype(float)
    ax.fill_between(x, p25_obs, p75_obs, color="#0984e3",
                    alpha=0.12, label="IQR obs (P25–P75)", zorder=1)

    # IQR modelo
    p25_mod = df_c["P25_mod"].values.astype(float)
    p75_mod = df_c["P75_mod"].values.astype(float)
    ax.fill_between(x, p25_mod, p75_mod, color="#fdcb6e",
                    alpha=0.12, label="IQR mod (P25–P75)", zorder=1)

    # Umbrales NOM-172
    for cat_key, cat_vals in NOM172.items():
        u_nom = cat_vals.get(cont)
        if u_nom is None: continue
        col_nom = "#FF7E00" if cat_key=="mala" else "#CC0000"
        ax.axhline(u_nom, color=col_nom, lw=1.2, ls="--", alpha=0.80,
                   label=f"NOM-172 {cat_key.replace('_',' ')} ({u_nom})")

    ax.set_xticks(x); ax.set_xticklabels(etqs, fontsize=8)
    ax.set_xlabel("Mes", fontsize=9)
    ax.set_ylabel(f"Concentración  [{meta['unidad']}]", fontsize=9)
    ax.set_title(
        f"Serie de percentiles — {meta['nombre']}  |  {ciudad}  |  {hor_label}\n"
        f"Línea continua = observado  ·  Línea punteada = modelo",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=7, ncol=3, loc="upper left",
              bbox_to_anchor=(0, -0.14), framealpha=0.90)
    ax.grid(True, color="#DDD", lw=0.5)
    ax.tick_params(labelsize=8)
    fig.text(0.01, 0.01,
             "WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | ICAyCC, UNAM",
             fontsize=6.5, style="italic", color="gray")

    plt.tight_layout(rect=[0, 0.10, 1, 1])
    nombre = (f"percentiles_{cont}_serie_{ciudad}_"
              f"{hor_label.replace(' ','')}.png")
    ruta   = os.path.join(dir_salida, nombre)
    os.makedirs(dir_salida, exist_ok=True)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# PROCESAMIENTO PRINCIPAL POR CONTAMINANTE
# ──────────────────────────────────────────────────────────────────────────────

def procesar_contaminante(datos: pd.DataFrame, cont: str,
                           cols_mod: List[str], pct_ref: int,
                           dir_salida: str,
                           serie_ciudades: Optional[List[str]]) -> pd.DataFrame:
    """
    Calcula la tabla de percentiles y genera todas las gráficas para un
    contaminante dado.
    """
    print(f"\n── {cont} ──────────────────────────────────────────")
    df_tabla = construir_tabla_percentiles(datos, cont, cols_mod)
    if df_tabla.empty:
        print(f"[WARN] {cont}: tabla vacía.")
        return pd.DataFrame()

    # CSV por contaminante
    ruta_csv = os.path.join(dir_salida, f"tabla_percentiles_{cont}.csv")
    os.makedirs(dir_salida, exist_ok=True)
    df_tabla.to_csv(ruta_csv, index=False)
    print(f"[OK]   CSV: {ruta_csv}")

    for col_mod in cols_mod:
        hor_label = HORIZONTES_LABELS[col_mod]
        print(f"  → {hor_label}")

        # Heatmap de percentil observado
        graficar_heatmap_obs(df_tabla, cont, pct_ref, col_mod, dir_salida)

        # Heatmap de sesgo percentílico
        graficar_heatmap_sesgo(df_tabla, cont, pct_ref, col_mod, dir_salida)

        # Diagrama Q-Q por temporada
        graficar_qq_temporada(datos, cont, col_mod, dir_salida)

        # Serie temporal por ciudad (solo ciudades con suficientes datos)
        ciudades_disponibles = sorted(
            df_tabla[df_tabla["horizonte"]==hor_label]["ciudad"].unique(),
            key=lambda c: CIUDADES_DOMINIO.index(c) if c in CIUDADES_DOMINIO else 99,
        )
        ciudades_serie = serie_ciudades or ciudades_disponibles[:3]
        for ciudad in ciudades_serie:
            if ciudad in ciudades_disponibles:
                graficar_serie_percentiles(
                    df_tabla, cont, col_mod, ciudad, dir_salida
                )

    return df_tabla


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="percentiles_obs.py",
        description=(
            "Calcula percentiles empíricos de las observaciones SINAICA y los "
            "compara contra el modelo WRF-Chem. Genera heatmaps de percentil "
            "observado, sesgo percentílico, diagramas Q-Q por temporada y "
            "series temporales de P25/P50/P75/P90."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--entrada",  "-i", default="combinado/ajustados",
                   help="Directorio con eval_*.csv (default: combinado/ajustados)")
    p.add_argument("--salida",   "-o", default="resultados_percentiles",
                   help="Directorio de salida (default: resultados_percentiles)")
    p.add_argument("--cont",     "-C", nargs="+", default=None,
                   choices=list(META_CONT.keys()),
                   help="Contaminante(s) a procesar (default: todos)")
    p.add_argument("--horizonte","-H", default="todos",
                   help="Horizonte: 24h, 48h, 72h o 'todos' (default: todos)")
    p.add_argument("--percentil","-p", type=int, default=90,
                   choices=PCTS,
                   help="Percentil de referencia para heatmaps (default: 90)")
    p.add_argument("--ciudades", "-c", nargs="+", default=None,
                   help=f"Ciudades a incluir. Catálogo: {', '.join(CIUDADES_DOMINIO)}")
    p.add_argument("--serie-ciudades", nargs="+", default=None,
                   help="Ciudades para las que generar la serie temporal "
                        "(default: las primeras 3 con datos)")
    p.add_argument("--inicio",   default=None, metavar="YYYY-MM",
                   help="Mes de inicio del período")
    p.add_argument("--fin",      default=None, metavar="YYYY-MM",
                   help="Mes de fin del período")
    p.add_argument("--umbral-pm25", type=float, default=UMBRAL_PM25,
                   help=f"Techo absoluto PM2.5 µg/m³ (default: {UMBRAL_PM25})")
    p.add_argument("--iqr-factor",  type=float, default=IQR_FACTOR,
                   help=f"Factor IQR outliers PM2.5 (default: {IQR_FACTOR})")
    p.add_argument("--dpi",      type=int,   default=DPI,
                   help=f"Resolución PNG (default: {DPI})")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# PUNTO DE ENTRADA
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    global IQR_FACTOR, UMBRAL_PM25, DPI
    IQR_FACTOR  = args.iqr_factor
    UMBRAL_PM25 = args.umbral_pm25
    DPI         = args.dpi

    # Horizontes a procesar
    if args.horizonte.lower() in ("todos", "all"):
        cols_mod = list(HORIZONTES_LABELS.keys())
    else:
        col = HORIZONTES_MAP.get(args.horizonte.lower().replace("+", ""))
        if col is None:
            sys.exit(f"[ERROR] Horizonte no reconocido: '{args.horizonte}'")
        cols_mod = [col]

    conts           = args.cont or list(META_CONT.keys())
    ciudades_filtro = parsear_lista_ciudades(args.ciudades)
    serie_ciudades  = parsear_lista_ciudades(args.serie_ciudades)

    print("=" * 66)
    print("  Análisis de Percentiles — WRF-Chem / SINAICA")
    print("=" * 66)
    print(f"  Entrada          : {args.entrada}")
    print(f"  Salida           : {args.salida}")
    print(f"  Contaminantes    : {conts}")
    print(f"  Horizontes       : {[HORIZONTES_LABELS[c] for c in cols_mod]}")
    print(f"  Percentil ref.   : P{args.percentil}")
    print(f"  Ciudades filtro  : {ciudades_filtro or 'todas'}")
    print(f"  Ciudades serie   : {serie_ciudades or 'primeras 3 con datos'}")
    print(f"  Período          : {args.inicio or 'inicio'} – {args.fin or 'fin'}")
    print(f"  QC PM2.5         : techo={UMBRAL_PM25} µg/m³  IQR·k={IQR_FACTOR}")
    print("=" * 66)

    datos = leer_datos(args.entrada, ciudades_filtro,
                       args.inicio, args.fin, None)

    tablas = []
    for cont in conts:
        if cont not in datos["contaminante"].unique():
            print(f"[WARN] {cont}: sin datos — omitido.")
            continue
        df_t = procesar_contaminante(
            datos, cont, cols_mod, args.percentil,
            args.salida, serie_ciudades,
        )
        if not df_t.empty:
            tablas.append(df_t)

    # CSV global con todos los contaminantes
    if tablas:
        df_global = pd.concat(tablas, ignore_index=True)
        ruta_global = os.path.join(args.salida, "tabla_percentiles_todos.csv")
        df_global.to_csv(ruta_global, index=False)
        print(f"\n[OK]  CSV global: {ruta_global}  ({len(df_global)} filas)")

    print(f"\n[DONE] Resultados guardados en '{args.salida}/'")


if __name__ == "__main__":
    main()
