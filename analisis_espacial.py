#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisis_espacial_percentil.py
================================
Evaluación dicotómica del pronóstico WRF-Chem usando umbrales adaptativos
locales derivados de los percentiles empíricos de las observaciones, en lugar
del umbral normativo nacional fijo (NOM-172-SEMARNAT-2023).

Motivación científica
---------------------
Un umbral fijo nacional (p.ej. O3 = 135 ppb) no refleja la climatología local
de cada ciudad. El P90 de O3 en Pachuca (≈145 ppb) es muy diferente al de
Cuernavaca (≈154 ppb). Evaluar el modelo con el mismo umbral en ambas ciudades
mezcla distribuciones distintas y puede distorsionar el diagnóstico:

  - Una ciudad con concentraciones basales ALTAS (Tula) tendrá muchos
    eventos verdaderos incluso si el modelo es mediocre.
  - Una ciudad con concentraciones basales BAJAS tendrá pocos eventos,
    lo que infla artificialmente el PC y el BIAS.

El análisis por percentiles locales responde la pregunta científicamente
más robusta:
    ¿El modelo detecta los días de concentración EXTREMA para CADA CIUDAD,
    independientemente del nivel absoluto de contaminación?

Metodología
-----------
  1. Para cada ciudad × mes × contaminante se calculan los percentiles
     empíricos de las OBSERVACIONES: P75, P90, P95 (configurables).
  2. Esos percentiles se usan como umbrales de evento en el cálculo
     dicotómico (POD, FAR, CSI, TSS, BIAS).
  3. Los resultados se comparan directamente contra la evaluación con
     umbral fijo NOM-172 para cuantificar la diferencia de diagnóstico.

Tres modos de umbral adaptativo
---------------------------------
  --modo ciudad_mes : umbral = Pxx de obs de esa ciudad en ese mes (más local)
  --modo ciudad     : umbral = Pxx de obs de esa ciudad en todo el período
  --modo dominio    : umbral = Pxx de obs de todas las ciudades (más global)

Salidas
-------
  espacial_<CONT>_<HOR>_<MODO>_P<N>.png
      Mapa de calor (ciudad × mes) de las cuatro métricas principales
      (POD, FAR, CSI, TSS) con el umbral adaptativo aplicado.

  espacial_<CONT>_<HOR>_comparacion_P<N>.png
      Panel comparativo: métricas con umbral adaptativo (izq.) vs
      umbral fijo NOM-172 (der.) para ver la diferencia de diagnóstico.

  espacial_<CONT>_<HOR>_umbral_adaptativo_P<N>.png
      Heatmap del valor del umbral adaptativo aplicado en cada celda
      (ciudad × mes), junto con la referencia NOM-172.

  tabla_espacial_percentil_<CONT>.csv
      Tabla completa: ciudad × mes × horizonte × percentil con
      umbral adaptativo, métricas dicotómicas y comparación vs NOM-172.

Uso
---
    # Todos los contaminantes, umbral P90 local por ciudad×mes
    python3 analisis_espacial_percentil.py

    # Solo O3, umbral P75 y P90, horizonte +24h
    python3 analisis_espacial_percentil.py --cont O3 --percentiles 75 90 --horizonte 24h

    # Modo ciudad (umbral estable por ciudad, sin variación mensual)
    python3 analisis_espacial_percentil.py --modo ciudad --cont O3

    # Comparar con NOM-172 en el mismo panel
    python3 analisis_espacial_percentil.py --cont O3 --comparar-nom172

    python3 analisis_espacial_percentil.py --help

Dependencias: pip install pandas numpy matplotlib scipy
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
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────

NOM172 = {
    "mala":     {"O3": 135.0, "PM10": 132.0, "PM25":  79.0, "SO2": 185.0},
    "muy_mala": {"O3": 175.0, "PM10": 213.0, "PM25": 130.0, "SO2": 304.0},
}

META_CONT = {
    "O3":   {"nombre": "Ozono (O₃)",              "unidad": "ppbv"},
    "PM10": {"nombre": "PM10",                    "unidad": "µg/m³"},
    "PM25": {"nombre": "PM2.5",                   "unidad": "µg/m³"},
    "SO2":  {"nombre": "Dióxido de azufre (SO₂)", "unidad": "ppbv"},
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
    "24h": "mod_dia1", "+24h": "mod_dia1",
    "48h": "mod_dia2", "+48h": "mod_dia2",
    "72h": "mod_dia3", "+72h": "mod_dia3",
    "todos": None,
}
HORIZONTES_LABELS = {
    "mod_dia1": "+24 h", "mod_dia2": "+48 h", "mod_dia3": "+72 h",
}

TEMPORADAS = {
    "Secas frías":     {"meses": {11, 12, 1, 2}, "color": "#4A90D9"},
    "Secas calientes": {"meses": {3, 4, 5},       "color": "#F5A623"},
    "Lluvias":         {"meses": {6, 7, 8, 9, 10},"color": "#2d8a2d"},
}

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
    1:"Ene", 2:"Feb",  3:"Mar",  4:"Abr",
    5:"May", 6:"Jun",  7:"Jul",  8:"Ago",
    9:"Sep", 10:"Oct", 11:"Nov", 12:"Dic",
}

# Métricas a mostrar en el heatmap principal
METRICAS_HEATMAP = ["POD", "FAR", "CSI", "TSS"]

# Semáforos de color por métrica (verde=bueno)
SEMAFORO_CFG = {
    "POD": {"cmap": "RdYlGn",   "vmin": 0.0, "vmax": 1.0, "ideal": 1.0},
    "FAR": {"cmap": "RdYlGn_r", "vmin": 0.0, "vmax": 1.0, "ideal": 0.0},
    "CSI": {"cmap": "RdYlGn",   "vmin": 0.0, "vmax": 1.0, "ideal": 1.0},
    "TSS": {"cmap": "RdYlGn",   "vmin":-1.0, "vmax": 1.0, "ideal": 1.0},
    "BIAS":{"cmap": "RdYlGn",   "vmin": 0.0, "vmax": 2.0, "ideal": 1.0},
}

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
    anio, mm = mes_str.split("-")
    return f"{NOMBRE_MES[int(mm)]}\n{anio[2:]}"


# ──────────────────────────────────────────────────────────────────────────────
# CONTROL DE CALIDAD
# ──────────────────────────────────────────────────────────────────────────────

def limpiar_serie(obs: np.ndarray, mod: np.ndarray,
                  cont: str) -> Tuple[np.ndarray, np.ndarray]:
    """NaN → límites físicos → filtro IQR (PM2.5)."""
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
    datos["fecha"]     = pd.to_datetime(datos["fecha"])
    datos["mes"]       = datos["fecha"].dt.month
    datos["mes_str"]   = datos["fecha"].dt.to_period("M").astype(str)
    datos["temporada"] = datos["mes"].map(mes_a_temporada)
    if inicio:
        datos = datos[datos["fecha"] >= pd.Timestamp(inicio+"-01")]
    if fin:
        datos = datos[datos["fecha"] <= pd.Timestamp(fin+"-01")+pd.offsets.MonthEnd(0)]
    print(f"[INFO] Registros: {len(datos):,}  "
          f"({datos['fecha'].dt.date.nunique()} días, "
          f"{datos['ciudad'].nunique()} ciudades)")
    return datos


# ──────────────────────────────────────────────────────────────────────────────
# CÁLCULO DEL UMBRAL ADAPTATIVO
# ──────────────────────────────────────────────────────────────────────────────

def calcular_umbral_adaptativo(
    datos_cont: pd.DataFrame,
    cont: str,
    pct: int,
    modo: str,
) -> Dict[Tuple, float]:
    """
    Calcula el umbral adaptativo para cada combinación (ciudad, mes_str)
    según el modo elegido:

        ciudad_mes : Pxx de obs de esa ciudad en ese mes   ← más local
        ciudad     : Pxx de obs de esa ciudad en todo el período
        dominio    : Pxx de obs de todas las ciudades juntas

    Retorna dict {(ciudad, mes_str): umbral}.
    """
    umbrales: Dict[Tuple, float] = {}
    ciudades  = datos_cont["ciudad"].unique()
    meses     = sorted(datos_cont["mes_str"].unique())

    if modo == "dominio":
        # Un único umbral para todo el dominio
        obs_all = np.asarray(datos_cont["max_obs"].values, float)
        obs_all = obs_all[np.isfinite(obs_all) & (obs_all > 0)]
        if len(obs_all) < MIN_DIAS:
            return umbrales
        u_dominio = float(np.percentile(obs_all, pct))
        for ciudad in ciudades:
            for mes_str in meses:
                umbrales[(ciudad, mes_str)] = u_dominio

    elif modo == "ciudad":
        # Un umbral por ciudad (todo el período)
        for ciudad in ciudades:
            dfc = datos_cont[datos_cont["ciudad"] == ciudad]
            obs = np.asarray(dfc["max_obs"].values, float)
            obs = obs[np.isfinite(obs) & (obs > 0)]
            if len(obs) < MIN_DIAS:
                continue
            u_ciudad = float(np.percentile(obs, pct))
            for mes_str in meses:
                umbrales[(ciudad, mes_str)] = u_ciudad

    else:   # ciudad_mes (más local)
        for ciudad in ciudades:
            dfc = datos_cont[datos_cont["ciudad"] == ciudad]
            for mes_str in meses:
                dfm = dfc[dfc["mes_str"] == mes_str]
                obs = np.asarray(dfm["max_obs"].values, float)
                obs = obs[np.isfinite(obs) & (obs > 0)]
                if len(obs) < MIN_DIAS:
                    continue
                umbrales[(ciudad, mes_str)] = float(np.percentile(obs, pct))

    return umbrales


# ──────────────────────────────────────────────────────────────────────────────
# CÁLCULO DE MÉTRICAS DICOTÓMICAS
# ──────────────────────────────────────────────────────────────────────────────

def contingencia(obs: np.ndarray, mod: np.ndarray,
                 umbral: float) -> Optional[Dict]:
    """Tabla de contingencia 2×2 y métricas para un umbral dado."""
    if len(obs) < MIN_DIAS:
        return None
    obs_ev = obs >= umbral
    mod_ev = mod >= umbral
    H = int(np.sum( obs_ev &  mod_ev))
    M = int(np.sum( obs_ev & ~mod_ev))
    F = int(np.sum(~obs_ev &  mod_ev))
    C = int(np.sum(~obs_ev & ~mod_ev))
    N = H + M + F + C
    N_ev = H + M

    def safe(num, den): return num / den if den > 0 else np.nan

    if N_ev == 0:
        return None     # sin eventos observados → métricas indefinidas

    POD  = safe(H, H+M)
    FAR  = safe(F, H+F)
    CSI  = safe(H, H+M+F)
    POFD = safe(F, F+C)
    TSS  = (POD - POFD) if np.isfinite(POD) and np.isfinite(POFD) else np.nan
    BIAS = safe(H+F, H+M)

    def r(x): return round(float(x), 4) if np.isfinite(x) else np.nan

    return {"H": H, "M": M, "F": F, "C": C,
            "N": N, "N_ev": N_ev,
            "POD": r(POD), "FAR": r(FAR), "CSI": r(CSI),
            "TSS": r(TSS), "BIAS": r(BIAS)}


# ──────────────────────────────────────────────────────────────────────────────
# CONSTRUCCIÓN DE LA TABLA PRINCIPAL
# ──────────────────────────────────────────────────────────────────────────────

def construir_tabla(
    datos: pd.DataFrame,
    cont: str,
    cols_mod: List[str],
    percentiles: List[int],
    modo: str,
) -> pd.DataFrame:
    """
    Construye el DataFrame completo con umbrales adaptativos y métricas
    para cada ciudad × mes × horizonte × percentil.

    Incluye también la comparación vs umbral fijo NOM-172 (categoría mala).
    """
    df_cont  = datos[datos["contaminante"] == cont].copy()
    ciudades = sorted(df_cont["ciudad"].unique(),
                      key=lambda c: CIUDADES_DOMINIO.index(c)
                      if c in CIUDADES_DOMINIO else 99)
    meses    = sorted(df_cont["mes_str"].unique())

    # Umbral fijo NOM-172 para comparación
    u_nom172 = NOM172["mala"].get(cont, np.nan)

    filas = []

    for pct in percentiles:
        print(f"    Calculando umbrales adaptativos P{pct} (modo={modo})...")
        umbrales_adapt = calcular_umbral_adaptativo(df_cont, cont, pct, modo)

        for col_mod in cols_mod:
            hor_label = HORIZONTES_LABELS[col_mod]

            for ciudad in ciudades:
                dfc = df_cont[df_cont["ciudad"] == ciudad]

                for mes_str in meses:
                    dfm = dfc[dfc["mes_str"] == mes_str]
                    if dfm.empty: continue

                    mes_num   = int(mes_str.split("-")[1])
                    temporada = mes_a_temporada(mes_num)

                    obs = np.asarray(dfm["max_obs"].values, float)
                    mod = np.asarray(dfm[col_mod].values, float)
                    obs, mod = limpiar_serie(obs, mod, cont)
                    if len(obs) < MIN_DIAS: continue

                    # Umbral adaptativo
                    u_adapt = umbrales_adapt.get((ciudad, mes_str))
                    if u_adapt is None: continue

                    st_adapt = contingencia(obs, mod, u_adapt)
                    if st_adapt is None: continue

                    # Métricas con umbral fijo NOM-172
                    st_nom = contingencia(obs, mod, u_nom172) \
                             if np.isfinite(u_nom172) else None

                    fila = {
                        "ciudad":         ciudad,
                        "mes_str":        mes_str,
                        "mes":            mes_num,
                        "temporada":      temporada,
                        "contaminante":   cont,
                        "horizonte":      hor_label,
                        "modo":           modo,
                        "percentil":      pct,
                        "umbral_adapt":   round(u_adapt, 2),
                        "umbral_nom172":  u_nom172,
                        # Diferencia absoluta entre umbral adaptativo y NOM-172
                        "dif_umbral":     round(u_adapt - u_nom172, 2)
                                          if np.isfinite(u_nom172) else np.nan,
                        # Métricas con umbral adaptativo
                        **{f"{k}_adapt": v for k, v in st_adapt.items()},
                    }

                    # Métricas con umbral NOM-172
                    if st_nom:
                        fila.update({f"{k}_nom172": v
                                     for k, v in st_nom.items()})
                    else:
                        for k in ["POD","FAR","CSI","TSS","BIAS",
                                  "H","M","F","C","N","N_ev"]:
                            fila[f"{k}_nom172"] = np.nan

                    # Delta: métrica_adapt − métrica_nom172
                    for metrica in ["POD","FAR","CSI","TSS"]:
                        va = fila.get(f"{metrica}_adapt", np.nan)
                        vn = fila.get(f"{metrica}_nom172", np.nan)
                        fila[f"delta_{metrica}"] = round(va-vn, 4) \
                            if (np.isfinite(va) and np.isfinite(vn)) else np.nan

                    filas.append(fila)

    df = pd.DataFrame(filas)
    print(f"[INFO] {cont}: {len(df)} filas en la tabla")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# HELPERS DE GRAFICACIÓN
# ──────────────────────────────────────────────────────────────────────────────

def _fondo_temporadas(ax, meses_str):
    for j, mes_str in enumerate(meses_str):
        mm = int(mes_str.split("-")[1])
        for nombre_t, cfg_t in TEMPORADAS.items():
            if mm in cfg_t["meses"]:
                ax.axvspan(j-0.5, j+0.5, color=cfg_t["color"],
                           alpha=0.15, zorder=0); break

def _formato_heatmap(ax, ciudades, meses_str, mat, cmap, vmin, vmax,
                     label_cb, fig, fmt=".2f", norm=None):
    """Dibuja el heatmap con anotaciones y barra de color."""
    if norm is None:
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(mat, cmap=cmap, aspect="auto",
                   norm=norm, interpolation="nearest", zorder=1)
    for i in range(len(ciudades)):
        for j in range(len(meses_str)):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, "N/D", ha="center", va="center",
                        fontsize=6.5, color="gray", zorder=3); continue
            rgba = im.cmap(norm(v))
            lum  = 0.299*rgba[0] + 0.587*rgba[1] + 0.114*rgba[2]
            col_txt = "white" if lum < 0.48 else "black"
            ax.text(j, i, f"{v:{fmt}}", ha="center", va="center",
                    fontsize=8, fontweight="bold", color=col_txt, zorder=3)
    ax.set_xticks(np.arange(-0.5, len(meses_str), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ciudades), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.9, zorder=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xticks(range(len(meses_str)))
    ax.set_xticklabels([etiqueta_mes(m) for m in meses_str], fontsize=7.5)
    ax.set_yticks(range(len(ciudades)))
    ax.set_yticklabels(ciudades, fontsize=8.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(label_cb, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return im

def _leyenda_temporadas(ax):
    parches = [mpatches.Patch(color=cfg["color"], alpha=0.65, label=n)
               for n, cfg in TEMPORADAS.items()]
    ax.legend(handles=parches, loc="upper left",
              bbox_to_anchor=(0, -0.14), ncol=3, fontsize=7, frameon=True)

def _nota_pie(fig):
    fig.text(0.01, 0.01,
             "WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | ICAyCC, UNAM",
             fontsize=6.5, style="italic", color="gray")


# ──────────────────────────────────────────────────────────────────────────────
# GRÁFICA 1 — HEATMAP 4 MÉTRICAS CON UMBRAL ADAPTATIVO
# ──────────────────────────────────────────────────────────────────────────────

def graficar_heatmap_metricas(df_t: pd.DataFrame, cont: str, col_mod: str,
                               pct: int, modo: str, dir_salida: str):
    """4 paneles (POD/FAR/CSI/TSS) con umbral adaptativo."""
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h = df_t[(df_t["horizonte"]==hor_label) & (df_t["percentil"]==pct)]
    if df_h.empty: return

    ciudades = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses    = sorted(df_h["mes_str"].unique())
    n_c, n_m = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    fig, axes = plt.subplots(2, 2, figsize=(max(10, n_m*0.9+3), n_c*0.85+3.5))

    for ax, metrica in zip(axes.flat, METRICAS_HEATMAP):
        cfg = SEMAFORO_CFG[metrica]
        mat = np.full((n_c, n_m), np.nan)
        for i, c in enumerate(ciudades):
            for j, m in enumerate(meses):
                row = df_h[(df_h["ciudad"]==c) & (df_h["mes_str"]==m)]
                if not row.empty:
                    mat[i, j] = row[f"{metrica}_adapt"].values[0]

        _fondo_temporadas(ax, meses)
        _formato_heatmap(ax, ciudades, meses, mat,
                         cfg["cmap"], cfg["vmin"], cfg["vmax"],
                         f"{metrica}  (P{pct} adaptativo)", fig)
        ax.set_title(f"{metrica}  —  ideal={cfg['ideal']:.0f}",
                     fontsize=10, fontweight="bold")
        _leyenda_temporadas(ax)

    fig.suptitle(
        f"Evaluación dicotómica con umbral P{pct} adaptativo ({modo})\n"
        f"{meta['nombre']}  |  Horizonte {hor_label}  |  "
        f"Umbral = percentil P{pct} de obs local",
        fontsize=11, fontweight="bold", y=1.02,
    )
    _nota_pie(fig)
    plt.tight_layout()
    nombre = (f"espacial_{cont}_{hor_label.replace(' ','')}"
              f"_{modo}_P{pct}.png")
    ruta = os.path.join(dir_salida, nombre)
    os.makedirs(dir_salida, exist_ok=True)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# GRÁFICA 2 — HEATMAP DEL UMBRAL ADAPTATIVO
# ──────────────────────────────────────────────────────────────────────────────

def graficar_heatmap_umbral(df_t: pd.DataFrame, cont: str, col_mod: str,
                             pct: int, modo: str, dir_salida: str):
    """
    Heatmap del valor numérico del umbral adaptativo aplicado en cada
    celda (ciudad × mes), con la referencia NOM-172 marcada en la
    barra de color.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h = df_t[(df_t["horizonte"]==hor_label) & (df_t["percentil"]==pct)]
    if df_h.empty: return

    ciudades = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses    = sorted(df_h["mes_str"].unique())
    n_c, n_m = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    mat_u  = np.full((n_c, n_m), np.nan)
    mat_dif = np.full((n_c, n_m), np.nan)
    u_nom  = df_h["umbral_nom172"].iloc[0] if "umbral_nom172" in df_h else np.nan

    for i, c in enumerate(ciudades):
        for j, m in enumerate(meses):
            row = df_h[(df_h["ciudad"]==c) & (df_h["mes_str"]==m)]
            if row.empty: continue
            mat_u[i, j]   = row["umbral_adapt"].values[0]
            mat_dif[i, j] = row["dif_umbral"].values[0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(14, n_m*1.6+4), n_c*0.85+3))

    # Panel izquierdo: valor absoluto del umbral adaptativo
    _fondo_temporadas(ax1, meses)
    vmin_u = np.nanmin(mat_u); vmax_u = np.nanmax(mat_u)
    im1 = _formato_heatmap(ax1, ciudades, meses, mat_u,
                            "YlOrRd", vmin_u, vmax_u,
                            f"Umbral P{pct} adaptativo  [{meta['unidad']}]",
                            fig, fmt=".1f")
    if np.isfinite(u_nom):
        norm1 = plt.Normalize(vmin_u, vmax_u)
        cb1   = fig.axes[-1]   # última colorbar añadida
        cb1.axhline((u_nom-vmin_u)/(vmax_u-vmin_u) if vmax_u>vmin_u else 0.5,
                    color="#FF7E00", lw=2.0, ls="--")
        cb1.text(1.6, (u_nom-vmin_u)/(vmax_u-vmin_u),
                 f"NOM-172\nmala={u_nom:.0f}",
                 va="center", fontsize=6.5, color="#FF7E00",
                 transform=cb1.transData)
    ax1.set_title(f"Umbral P{pct} adaptativo por ciudad × mes",
                  fontsize=10, fontweight="bold")
    _leyenda_temporadas(ax1)

    # Panel derecho: diferencia umbral_adapt − umbral_nom172
    _fondo_temporadas(ax2, meses)
    lim_dif = np.nanpercentile(np.abs(mat_dif[np.isfinite(mat_dif)]), 95) \
              if np.any(np.isfinite(mat_dif)) else 10
    lim_dif = max(lim_dif, 1.0)
    norm_div = TwoSlopeNorm(vcenter=0.0, vmin=-lim_dif, vmax=lim_dif)
    _formato_heatmap(ax2, ciudades, meses, mat_dif,
                     "RdBu", -lim_dif, lim_dif,
                     f"P{pct} adapt − NOM-172  [{meta['unidad']}]",
                     fig, fmt=".1f", norm=norm_div)
    ax2.set_title(
        f"Diferencia umbral adaptativo − NOM-172 mala\n"
        f"Azul = P{pct} > NOM (umbral más alto)  ·  "
        f"Rojo = P{pct} < NOM (umbral más bajo)",
        fontsize=9, fontweight="bold")
    _leyenda_temporadas(ax2)

    fig.suptitle(
        f"Umbral adaptativo P{pct} ({modo}) — {meta['nombre']}  |  "
        f"Horizonte {hor_label}",
        fontsize=11, fontweight="bold", y=1.02,
    )
    _nota_pie(fig)
    plt.tight_layout()
    nombre = (f"espacial_{cont}_{hor_label.replace(' ','')}"
              f"_umbral_P{pct}_{modo}.png")
    ruta = os.path.join(dir_salida, nombre)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# GRÁFICA 3 — COMPARACIÓN UMBRAL ADAPTATIVO vs NOM-172
# ──────────────────────────────────────────────────────────────────────────────

def graficar_comparacion(df_t: pd.DataFrame, cont: str, col_mod: str,
                          pct: int, modo: str, dir_salida: str):
    """
    Panel 2×4: cada fila es una métrica (POD, FAR, CSI, TSS).
    Columna izquierda = umbral adaptativo P<pct>.
    Columna derecha   = umbral fijo NOM-172.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h = df_t[(df_t["horizonte"]==hor_label) & (df_t["percentil"]==pct)]
    if df_h.empty: return

    ciudades = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses    = sorted(df_h["mes_str"].unique())
    n_c, n_m = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    n_met  = len(METRICAS_HEATMAP)
    fig, axes = plt.subplots(n_met, 2,
                              figsize=(max(14, n_m*1.6+4), n_c*n_met*0.7+2))
    if n_met == 1: axes = axes.reshape(1, 2)

    for row_idx, metrica in enumerate(METRICAS_HEATMAP):
        cfg  = SEMAFORO_CFG[metrica]
        ax_a = axes[row_idx, 0]   # umbral adaptativo
        ax_n = axes[row_idx, 1]   # NOM-172

        mat_a = np.full((n_c, n_m), np.nan)
        mat_n = np.full((n_c, n_m), np.nan)

        for i, c in enumerate(ciudades):
            for j, m in enumerate(meses):
                row = df_h[(df_h["ciudad"]==c) & (df_h["mes_str"]==m)]
                if row.empty: continue
                mat_a[i, j] = row[f"{metrica}_adapt"].values[0]
                mat_n[i, j] = row[f"{metrica}_nom172"].values[0]

        _fondo_temporadas(ax_a, meses)
        _fondo_temporadas(ax_n, meses)

        _formato_heatmap(ax_a, ciudades, meses, mat_a,
                         cfg["cmap"], cfg["vmin"], cfg["vmax"],
                         f"{metrica}  (P{pct} adaptativo)", fig)
        _formato_heatmap(ax_n, ciudades, meses, mat_n,
                         cfg["cmap"], cfg["vmin"], cfg["vmax"],
                         f"{metrica}  (NOM-172 mala)", fig)

        ax_a.set_title(f"{metrica}  —  P{pct} adaptativo ({modo})",
                       fontsize=9, fontweight="bold")
        ax_n.set_title(f"{metrica}  —  NOM-172 mala (umbral fijo)",
                       fontsize=9, fontweight="bold")

        if row_idx == n_met - 1:
            _leyenda_temporadas(ax_a)
            _leyenda_temporadas(ax_n)

    fig.suptitle(
        f"Comparación umbral adaptativo P{pct} vs NOM-172 — "
        f"{meta['nombre']}  |  Horizonte {hor_label}",
        fontsize=11, fontweight="bold", y=1.01,
    )
    _nota_pie(fig)
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.subplots_adjust(hspace=0.55, wspace=0.35)
    nombre = (f"espacial_{cont}_{hor_label.replace(' ','')}"
              f"_comparacion_P{pct}_{modo}.png")
    ruta = os.path.join(dir_salida, nombre)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# GRÁFICA 4 — DELTA (ADAPT − NOM172) POR MÉTRICA
# ──────────────────────────────────────────────────────────────────────────────

def graficar_delta(df_t: pd.DataFrame, cont: str, col_mod: str,
                   pct: int, modo: str, dir_salida: str):
    """
    Heatmap de la diferencia (métrica_adapt − métrica_nom172) para
    revelar dónde cambia el diagnóstico al usar el umbral adaptativo.
    Verde = el diagnóstico mejora; rojo = empeora.
    """
    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_h = df_t[(df_t["horizonte"]==hor_label) & (df_t["percentil"]==pct)]
    if df_h.empty: return

    ciudades = [c for c in CIUDADES_DOMINIO if c in df_h["ciudad"].unique()]
    meses    = sorted(df_h["mes_str"].unique())
    n_c, n_m = len(ciudades), len(meses)
    if n_c == 0 or n_m == 0: return

    # Métricas donde delta>0 es bueno: POD, CSI, TSS
    # Para FAR delta<0 es bueno → invertir escala
    fig, axes = plt.subplots(2, 2,
                              figsize=(max(10, n_m*0.9+3), n_c*0.85+3.5))

    for ax, metrica in zip(axes.flat, METRICAS_HEATMAP):
        mat = np.full((n_c, n_m), np.nan)
        for i, c in enumerate(ciudades):
            for j, m in enumerate(meses):
                row = df_h[(df_h["ciudad"]==c) & (df_h["mes_str"]==m)]
                if row.empty: continue
                mat[i, j] = row[f"delta_{metrica}"].values[0]

        lim = np.nanpercentile(np.abs(mat[np.isfinite(mat)]), 95) \
              if np.any(np.isfinite(mat)) else 0.1
        lim = max(lim, 0.01)

        # Para FAR: rojo=delta>0 (FAR sube=empeora); verde=delta<0
        cmap = "RdYlGn_r" if metrica == "FAR" else "RdYlGn"
        norm_div = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)

        _fondo_temporadas(ax, meses)
        _formato_heatmap(ax, ciudades, meses, mat, cmap, -lim, lim,
                         f"Δ{metrica}  (adapt − NOM-172)", fig,
                         fmt="+.3f", norm=norm_div)

        mejora_txt = "↑ = mejora" if metrica != "FAR" else "↓ = mejora"
        ax.set_title(f"Δ{metrica}  ({mejora_txt})",
                     fontsize=10, fontweight="bold")
        _leyenda_temporadas(ax)

    fig.suptitle(
        f"Delta diagnóstico: P{pct} adaptativo ({modo}) − NOM-172 mala\n"
        f"{meta['nombre']}  |  Horizonte {hor_label}  |  "
        f"Verde = diagnóstico mejora con umbral adaptativo",
        fontsize=11, fontweight="bold", y=1.02,
    )
    _nota_pie(fig)
    plt.tight_layout()
    nombre = (f"espacial_{cont}_{hor_label.replace(' ','')}"
              f"_delta_P{pct}_{modo}.png")
    ruta = os.path.join(dir_salida, nombre)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# PROCESAMIENTO POR CONTAMINANTE
# ──────────────────────────────────────────────────────────────────────────────

def procesar_cont(datos: pd.DataFrame, cont: str, cols_mod: List[str],
                  percentiles: List[int], modo: str,
                  comparar: bool, dir_salida: str) -> pd.DataFrame:

    print(f"\n── {cont} {'─'*50}")
    df_t = construir_tabla(datos, cont, cols_mod, percentiles, modo)
    if df_t.empty:
        print(f"[WARN] {cont}: tabla vacía.")
        return pd.DataFrame()

    ruta_csv = os.path.join(dir_salida, f"tabla_espacial_percentil_{cont}.csv")
    os.makedirs(dir_salida, exist_ok=True)
    df_t.to_csv(ruta_csv, index=False)
    print(f"[OK]   CSV: {ruta_csv}")

    for col_mod in cols_mod:
        for pct in percentiles:
            hor_label = HORIZONTES_LABELS[col_mod]
            print(f"  → {hor_label}  P{pct}")
            graficar_heatmap_metricas(df_t, cont, col_mod, pct, modo, dir_salida)
            graficar_heatmap_umbral  (df_t, cont, col_mod, pct, modo, dir_salida)
            graficar_delta           (df_t, cont, col_mod, pct, modo, dir_salida)
            if comparar:
                graficar_comparacion (df_t, cont, col_mod, pct, modo, dir_salida)

    return df_t


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="analisis_espacial_percentil.py",
        description=(
            "Evaluación dicotómica WRF-Chem con umbrales adaptativos locales "
            "derivados de percentiles empíricos de las observaciones SINAICA. "
            "Compara el diagnóstico con umbral adaptativo vs umbral fijo NOM-172."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--entrada",  "-i", default="combinado/ajustados")
    p.add_argument("--salida",   "-o", default="resultados_espacial_percentil")
    p.add_argument("--cont",     "-C", nargs="+", default=None,
                   choices=list(META_CONT.keys()))
    p.add_argument("--horizonte","-H", default="todos",
                   help="24h, 48h, 72h o 'todos' (default: todos)")
    p.add_argument("--percentiles","-p", nargs="+", type=int, default=[90],
                   help="Percentil(es) para el umbral adaptativo (default: 90)")
    p.add_argument("--modo",     "-m", default="ciudad_mes",
                   choices=["ciudad_mes","ciudad","dominio"],
                   help=(
                       "Ámbito del umbral adaptativo:\n"
                       "  ciudad_mes — P<N> de obs de esa ciudad en ese mes (default)\n"
                       "  ciudad     — P<N> de obs de esa ciudad en todo el período\n"
                       "  dominio    — P<N> de obs de todas las ciudades"
                   ))
    p.add_argument("--comparar-nom172", action="store_true",
                   help="Generar panel comparativo adapt vs NOM-172")
    p.add_argument("--ciudades", "-c", nargs="+", default=None)
    p.add_argument("--inicio",   default=None, metavar="YYYY-MM")
    p.add_argument("--fin",      default=None, metavar="YYYY-MM")
    p.add_argument("--umbral-pm25", type=float, default=UMBRAL_PM25)
    p.add_argument("--iqr-factor",  type=float, default=IQR_FACTOR)
    p.add_argument("--dpi",      type=int, default=DPI)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    global IQR_FACTOR, UMBRAL_PM25, DPI
    IQR_FACTOR  = args.iqr_factor
    UMBRAL_PM25 = args.umbral_pm25
    DPI         = args.dpi

    if args.horizonte.lower() in ("todos","all"):
        cols_mod = list(HORIZONTES_LABELS.keys())
    else:
        col = HORIZONTES_MAP.get(args.horizonte.lower().replace("+",""))
        if col is None:
            sys.exit(f"[ERROR] Horizonte no reconocido: '{args.horizonte}'")
        cols_mod = [col]

    conts           = args.cont or list(META_CONT.keys())
    ciudades_filtro = parsear_lista_ciudades(args.ciudades)
    percentiles     = sorted(set(args.percentiles))

    print("=" * 66)
    print("  Análisis Espacial por Percentil — WRF-Chem / SINAICA")
    print("=" * 66)
    print(f"  Entrada       : {args.entrada}")
    print(f"  Salida        : {args.salida}")
    print(f"  Contaminantes : {conts}")
    print(f"  Horizontes    : {[HORIZONTES_LABELS[c] for c in cols_mod]}")
    print(f"  Percentil(es) : {percentiles}")
    print(f"  Modo umbral   : {args.modo}")
    print(f"  Comparar NOM  : {args.comparar_nom172}")
    print(f"  Ciudades      : {ciudades_filtro or 'todas'}")
    print(f"  QC PM2.5      : techo={UMBRAL_PM25}  IQR*k={IQR_FACTOR}")
    print("=" * 66)

    datos = leer_datos(args.entrada, ciudades_filtro,
                       args.inicio, args.fin, None)

    tablas = []
    for cont in conts:
        if cont not in datos["contaminante"].unique():
            print(f"[WARN] {cont}: sin datos — omitido.")
            continue
        df_t = procesar_cont(datos, cont, cols_mod, percentiles,
                             args.modo, args.comparar_nom172, args.salida)
        if not df_t.empty:
            tablas.append(df_t)

    if tablas:
        df_g = pd.concat(tablas, ignore_index=True)
        ruta_g = os.path.join(args.salida, "tabla_espacial_todos.csv")
        df_g.to_csv(ruta_g, index=False)
        print(f"\n[OK]  CSV global: {ruta_g}  ({len(df_g)} filas)")

    print(f"\n[DONE] Resultados en '{args.salida}/'")


if __name__ == "__main__":
    main()

