#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curva_habilidad.py
==================
Analiza la sensibilidad del pronóstico WRF-Chem a distintos umbrales de
concentración mediante el barrido continuo de niveles de alerta.

Para cada contaminante y ciudad se recorre un rango de umbrales y se calculan
las métricas dicotómicas POD, SR (= 1 − FAR), CSI, TSS y BIAS en cada nivel.
El resultado permite identificar:

    1. El umbral donde el modelo maximiza su habilidad (CSI_max)
    2. Si el umbral normativo de la NOM-172 cae en la zona óptima
    3. El comportamiento estructural del modelo: sesgo bajo vs alto umbral
    4. Diferencias estacionales (secas vs lluvias) en la habilidad

Salidas por contaminante
------------------------
    curva_habilidad_<CONT>_<HOR>_<AGRUP>.png
        Panel con 4 subgráficas:
          (a) POD y SR vs umbral  — con banda de incertidumbre binomial 90 %
          (b) CSI y TSS vs umbral — con marcador del CSI_max
          (c) BIAS de frecuencia vs umbral
          (d) Trayectoria de habilidad en el espacio Roebber (multi-umbral)

    tabla_barrido_umbrales.csv
        Tabla completa: umbral × agrupación × horizonte → todas las métricas

Fundamento estadístico
----------------------
    La "trayectoria de habilidad" en el espacio Roebber (SR, POD) al variar el
    umbral u ∈ [u_min, u_max] es una curva paramétrica coloreada por nivel de
    umbral (plasma_r). Su forma revela el comportamiento estructural del modelo:

    • Curva cercana a la diagonal b=1  → sin sesgo de frecuencia independiente
      del umbral elegido.
    • Curva que cruza la diagonal      → el modelo pasa de subestimar a
      sobreestimar la frecuencia de eventos al cambiar el umbral.
    • CSI_max distante del umbral NOM  → el nivel de alerta óptimo del modelo
      difiere del normativo (oportunidad de recalibración).

Umbrales normativos marcados — NOM-172-SEMARNAT-2023
-----------------------------------------------------
    Categoría Mala     : O3=135ppb  PM10=132µg/m³  PM25=79µg/m³   SO2=185ppb
    Categoría Muy Mala : O3=175ppb  PM10=213µg/m³  PM25=130µg/m³  SO2=304ppb

Uso
---
    # Todos los contaminantes, agrupado por ciudad
    python3 curva_habilidad.py

    # Solo O3, horizonte +24h, agrupado por temporada
    python3 curva_habilidad.py --cont O3 --horizonte 24h --agrupar temporada

    # Agrupado por horizonte (todos en un diagrama)
    python3 curva_habilidad.py --cont O3 --agrupar horizonte

    # Paso fino, mínimo de eventos por punto
    python3 curva_habilidad.py --paso 2 --min-eventos 5

    # Filtro de ciudades y período
    python3 curva_habilidad.py --ciudades CDMX Tula Pachuca --inicio 2026-01

    python3 curva_habilidad.py --help

Dependencias
------------
    pip install pandas numpy matplotlib scipy

Autor  : Pipeline ddsinaica / WRF-Chem — ICAyCC, UNAM
Versión: 1.0.0 (2026-07)

Referencias
-----------
    Roebber, P.J., 2009: Wea. Forecasting, 24, 749-755.
    Mason & Graham, 2002: Q.J.R. Meteorol. Soc., 128, 2145-2166.
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
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN GLOBAL
# ──────────────────────────────────────────────────────────────────────────────

NOM172 = {
    "mala":     {"O3": 135.0, "PM10": 132.0, "PM25":  79.0, "SO2": 185.0},
    "muy_mala": {"O3": 175.0, "PM10": 213.0, "PM25": 130.0, "SO2": 304.0},
}

RANGOS_DEFAULT = {
    "O3":   (40.0,  180.0, 5.0),
    "PM10": (10.0,  220.0, 5.0),
    "PM25": ( 5.0,  140.0, 5.0),
    "SO2":  (10.0,  320.0, 5.0),
}

META_CONT = {
    "O3":   {"nombre": "Ozono (O3)",              "unidad": "ppbv",  "color": "#1f77b4"},
    "PM10": {"nombre": "PM10",                    "unidad": "ug/m3", "color": "#d62728"},
    "PM25": {"nombre": "PM2.5",                   "unidad": "ug/m3", "color": "#ff7f0e"},
    "SO2":  {"nombre": "Dioxido de azufre (SO2)", "unidad": "ppbv",  "color": "#9467bd"},
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

COLORES_CIUDAD = {
    "CDMX": "#e6194b", "Cuernavaca": "#f58231", "Pachuca": "#bfef45",
    "Puebla": "#3cb44b", "SJdelRio": "#42d4f4", "Tlaxcala": "#4363d8",
    "Toluca": "#911eb4", "Tula": "#a9a9a9",
}
COLORES_HOR = {
    "mod_dia1": "#1f77b4", "mod_dia2": "#ff7f0e", "mod_dia3": "#d62728",
}

TEMPORADAS = {
    "Secas frias":    {"meses": {11,12,1,2}, "color": "#4A90D9"},
    "Secas calientes":{"meses": {3,4,5},     "color": "#F5A623"},
    "Lluvias":        {"meses": {6,7,8,9,10},"color": "#2d8a2d"},
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
MIN_EVENTOS = 3
DPI         = 150


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
                 f"        Validas: {CIUDADES_DOMINIO}")
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

def _smooth(arr: np.ndarray, window: int = 7) -> np.ndarray:
    if len(arr) < window or window < 3: return arr
    try:
        wl = min(window, len(arr))
        wl = wl if wl % 2 == 1 else wl - 1
        return savgol_filter(arr, window_length=wl, polyorder=2)
    except Exception:
        return arr


# ──────────────────────────────────────────────────────────────────────────────
# CONTROL DE CALIDAD
# ──────────────────────────────────────────────────────────────────────────────

def limpiar_serie(obs: np.ndarray, mod: np.ndarray,
                  cont: str) -> Tuple[np.ndarray, np.ndarray]:
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
# CALCULO DE METRICAS
# ──────────────────────────────────────────────────────────────────────────────

def contingencia_escalar(obs: np.ndarray, mod: np.ndarray,
                         umbral: float) -> Dict:
    obs_ev = obs >= umbral
    mod_ev = mod >= umbral
    H = int(np.sum( obs_ev &  mod_ev))
    M = int(np.sum( obs_ev & ~mod_ev))
    F = int(np.sum(~obs_ev &  mod_ev))
    C = int(np.sum(~obs_ev & ~mod_ev))
    N        = H + M + F + C
    N_ev_obs = H + M
    N_ev_mod = H + F

    def safe(num, den): return num / den if den > 0 else np.nan

    nan_row = {"POD": np.nan, "SR": np.nan, "FAR": np.nan,
               "CSI": np.nan, "TSS": np.nan, "BIAS": np.nan,
               "H": H, "M": M, "F": F, "C": C,
               "N": N, "N_ev_obs": N_ev_obs, "N_ev_mod": N_ev_mod}

    if N_ev_obs < MIN_EVENTOS:
        return nan_row

    POD  = safe(H, H+M)
    FAR  = safe(F, H+F)
    SR   = (1.0 - FAR) if np.isfinite(FAR) else np.nan
    CSI  = safe(H, H+M+F)
    POFD = safe(F, F+C)
    TSS  = (POD - POFD) if (np.isfinite(POD) and np.isfinite(POFD)) else np.nan
    BIAS = safe(H+F, H+M)

    return {"POD": POD, "SR": SR, "FAR": FAR, "CSI": CSI,
            "TSS": TSS, "BIAS": BIAS,
            "H": H, "M": M, "F": F, "C": C,
            "N": N, "N_ev_obs": N_ev_obs, "N_ev_mod": N_ev_mod}


def barrido_umbrales(obs: np.ndarray, mod: np.ndarray,
                     u_min: float, u_max: float, paso: float) -> pd.DataFrame:
    umbrales = np.arange(u_min, u_max + paso*0.5, paso)
    filas = []
    for u in umbrales:
        row = contingencia_escalar(obs, mod, u)
        row["umbral"]   = round(float(u), 4)
        row["frec_obs"] = row["N_ev_obs"] / row["N"] if row["N"] > 0 else np.nan
        row["frec_mod"] = row["N_ev_mod"] / row["N"] if row["N"] > 0 else np.nan
        filas.append(row)
    return pd.DataFrame(filas)


def calcular_csi_max(df: pd.DataFrame) -> Tuple[float, float]:
    df_v = df.dropna(subset=["CSI"])
    if df_v.empty: return np.nan, np.nan
    idx = df_v["CSI"].idxmax()
    return float(df_v.loc[idx, "umbral"]), float(df_v.loc[idx, "CSI"])


def rango_automatico(obs_global: np.ndarray, cont: str,
                     paso: float) -> Tuple[float, float]:
    d_min, d_max, _ = RANGOS_DEFAULT.get(cont, (0, 300, 5))
    if len(obs_global) < 10: return d_min, d_max
    p5  = float(np.percentile(obs_global,  5))
    p99 = float(np.percentile(obs_global, 99))
    margen = (p99 - p5) * 0.10
    u_min = max(d_min, p5  - margen)
    u_max = min(d_max, p99 + margen)
    return float(np.floor(u_min/paso)*paso), float(np.ceil(u_max/paso)*paso)


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
        sys.exit("[ERROR] Ningun archivo coincidio con los filtros.")
    datos = pd.concat(partes, ignore_index=True).rename(columns={"Fecha":"fecha"})
    datos["fecha"] = pd.to_datetime(datos["fecha"])
    if inicio:
        datos = datos[datos["fecha"] >= pd.Timestamp(inicio + "-01")]
    if fin:
        datos = datos[datos["fecha"] <= pd.Timestamp(fin+"-01")+pd.offsets.MonthEnd(0)]
    print(f"[INFO] Registros cargados: {len(datos):,} "
          f"({datos['fecha'].dt.date.nunique()} dias unicos)")
    return datos


# ──────────────────────────────────────────────────────────────────────────────
# FONDO DEL DIAGRAMA DE ROEBBER
# ──────────────────────────────────────────────────────────────────────────────

def _fondo_roebber(ax: plt.Axes):
    sr_g  = np.linspace(0.001, 1.0, 400)
    pod_g = np.linspace(0.001, 1.0, 400)
    SR2, POD2 = np.meshgrid(sr_g, pod_g)
    CSI2 = np.clip(1.0/(1.0/SR2 + 1.0/POD2 - 1.0), 0, 1)
    ax.contourf(SR2, POD2, CSI2, levels=np.linspace(0,1,21),
                cmap="YlGn", alpha=0.15, zorder=0)
    cs = ax.contour(SR2, POD2, CSI2,
                    levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],
                    colors="#2d6a2d", linewidths=0.7, alpha=0.6, zorder=1)
    ax.clabel(cs, fmt="%.1f", fontsize=6.5, inline=True, colors="#2d6a2d")
    for bv in [0.5, 0.75, 1.0, 1.5, 2.0]:
        sr_l  = np.linspace(0, 1, 200)
        pod_l = bv * sr_l
        mask  = (pod_l >= 0) & (pod_l <= 1)
        lw  = 1.5 if bv == 1.0 else 0.7
        col = "#333" if bv == 1.0 else "#999"
        ax.plot(sr_l[mask], pod_l[mask], ls="--", lw=lw, color=col,
                alpha=0.6, zorder=1)
        ie = np.where(mask)[0]
        if len(ie):
            ax.text(sr_l[ie[-1]]+0.01, pod_l[ie[-1]],
                    f"b={bv:.2g}", fontsize=6, color=col, va="bottom")
    ax.plot(1, 1, "*", ms=14, color="gold", markeredgecolor="#333",
            markeredgewidth=0.8, zorder=8, label="Pronostico perfecto")


# ──────────────────────────────────────────────────────────────────────────────
# GRAFICACION PRINCIPAL
# ──────────────────────────────────────────────────────────────────────────────

def graficar_curva_habilidad(
    resultados: Dict,
    cont: str,
    titulo: str,
    dir_salida: str,
    nombre_archivo: str,
    suavizar: bool = True,
    paleta: Optional[Dict[str, str]] = None,
):
    meta   = META_CONT[cont]
    unidad = meta["unidad"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_pod, ax_csi = axes[0,0], axes[0,1]
    ax_bias, ax_rob = axes[1,0], axes[1,1]

    # Configuracion basica de ejes
    for ax, ylabel, title in [
        (ax_pod,  "Valor",                     "(a) POD y SR vs umbral"),
        (ax_csi,  "Valor",                     "(b) CSI y TSS vs umbral"),
        (ax_bias, f"BIAS freq  (ideal = 1)",   "(c) BIAS de frecuencia vs umbral"),
    ]:
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel(f"Umbral  [{unidad}]", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.grid(True, color="#DDD", lw=0.5)
        ax.tick_params(labelsize=8)

    ax_pod.set_ylim(-0.05, 1.10)
    ax_csi.set_ylim(-0.18, 1.10)
    ax_pod.axhline(0, color="#ccc", lw=0.5, ls=":")
    ax_pod.axhline(1, color="#ccc", lw=0.5, ls=":")
    ax_csi.axhline(0, color="#ccc", lw=0.5, ls=":")
    ax_bias.axhline(1.0, color="#333", lw=1.6, ls="--", alpha=0.75,
                    label="BIAS = 1  (sin sesgo de frecuencia)", zorder=3)
    ax_bias.axhline(0, color="#ccc", lw=0.5, ls=":")

    # Marcar umbrales NOM-172 en los tres ejes lineales
    for cat_key, cat_vals in NOM172.items():
        u_nom = cat_vals.get(cont)
        if u_nom is None: continue
        col   = "#FF7E00" if cat_key == "mala" else "#CC0000"
        label = f"NOM-172 {cat_key.replace('_',' ')} ({u_nom} {unidad})"
        for ax in [ax_pod, ax_csi, ax_bias]:
            ax.axvline(u_nom, color=col, lw=1.3, ls="--",
                       alpha=0.80, zorder=4, label=label)

    # Diagrama de Roebber
    _fondo_roebber(ax_rob)
    ax_rob.set_title("(d) Trayectoria de habilidad multi-umbral (Roebber)",
                     fontsize=10, fontweight="bold")
    ax_rob.set_xlabel("SR  (1 - FAR = H / (H + F))", fontsize=9)
    ax_rob.set_ylabel("POD  (H / (H + M))", fontsize=9)
    ax_rob.set_xlim(-0.02, 1.05); ax_rob.set_ylim(-0.02, 1.05)
    ax_rob.set_aspect("equal")
    ax_rob.grid(True, color="#DDD", lw=0.5)
    ax_rob.tick_params(labelsize=8)

    cmap_tray = plt.colormaps["plasma_r"]

    for etiq, df_b in resultados.items():
        color = (paleta or {}).get(etiq, "#1f77b4")
        df_v  = df_b.dropna(subset=["POD","SR"])
        if df_v.empty: continue

        umbrales = df_v["umbral"].values
        pod  = df_v["POD"].values;  sr   = df_v["SR"].values
        csi  = df_v["CSI"].fillna(0).values
        tss  = df_v["TSS"].fillna(0).values
        bias = df_v["BIAS"].values
        n_tot= df_v["N"].values

        if suavizar and len(pod) >= 7:
            pod_s  = _smooth(pod, 7);  sr_s   = _smooth(sr,  7)
            csi_s  = _smooth(csi, 7);  tss_s  = _smooth(tss, 7)
            bias_s = _smooth(bias, 7)
        else:
            pod_s=pod; sr_s=sr; csi_s=csi; tss_s=tss; bias_s=bias

        u_opt, csi_opt = calcular_csi_max(df_v)

        # ─ (a) POD y SR ───────────────────────────────────────────────────────
        ax_pod.plot(umbrales, pod_s, color=color, lw=2.0,
                    ls="-",  label=f"{etiq} — POD")
        ax_pod.plot(umbrales, sr_s,  color=color, lw=2.0,
                    ls="--", alpha=0.75, label=f"{etiq} — SR")
        # Banda de incertidumbre binomial al 90 %
        se_pod = np.sqrt(np.clip(pod,0,1)*(1-np.clip(pod,0,1)) /
                         np.where(n_tot>0, n_tot, 1)) * 1.645
        ax_pod.fill_between(umbrales,
                            np.clip(pod_s - se_pod, 0, 1),
                            np.clip(pod_s + se_pod, 0, 1),
                            color=color, alpha=0.12)

        # ─ (b) CSI y TSS ──────────────────────────────────────────────────────
        ax_csi.plot(umbrales, csi_s, color=color, lw=2.0,
                    ls="-",  label=f"{etiq} — CSI")
        ax_csi.plot(umbrales, tss_s, color=color, lw=2.0,
                    ls=":",  alpha=0.80, label=f"{etiq} — TSS")
        # Marcar CSI_max
        if np.isfinite(u_opt) and np.isfinite(csi_opt):
            ax_csi.axvline(u_opt, color=color, lw=1.0, ls=":", alpha=0.7, zorder=5)
            ax_csi.annotate(
                f"CSI$_{{max}}$={csi_opt:.2f}\n@{u_opt:.0f}",
                xy=(u_opt, csi_opt),
                xytext=(u_opt + (umbrales[-1]-umbrales[0])*0.04,
                        max(0.0, csi_opt - 0.10)),
                fontsize=7, color=color,
                arrowprops=dict(arrowstyle="->", color=color, lw=0.8),
                zorder=6,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
            )

        # ─ (c) BIAS ───────────────────────────────────────────────────────────
        bias_plot = np.clip(bias_s, 0, 5)
        ax_bias.plot(umbrales, bias_plot, color=color, lw=2.0, label=etiq)
        ax_bias.set_ylim(0, min(5, np.nanpercentile(bias[np.isfinite(bias)], 97)*1.3)
                         if np.any(np.isfinite(bias)) else 5)

        # ─ (d) Trayectoria de Roebber ─────────────────────────────────────────
        if len(sr_s) > 1:
            norm_u = plt.Normalize(umbrales.min(), umbrales.max())
            for j in range(len(sr_s) - 1):
                c_seg = cmap_tray(norm_u(umbrales[j]))
                ax_rob.plot([sr_s[j], sr_s[j+1]],
                            [pod_s[j], pod_s[j+1]],
                            color=c_seg, lw=2.2, alpha=0.88, zorder=3)
            ax_rob.scatter(sr_s[0],  pod_s[0],  c=[[cmap_tray(0.05)]],
                           s=70, zorder=5, edgecolors="white", lw=0.8,
                           marker="o", label=f"{etiq} (u minimo)")
            ax_rob.scatter(sr_s[-1], pod_s[-1], c=[[cmap_tray(0.95)]],
                           s=70, zorder=5, edgecolors="white", lw=0.8,
                           marker="s", label=f"{etiq} (u maximo)")

            # Marcar umbrales NOM-172 en la trayectoria
            for cat_key, cat_vals in NOM172.items():
                u_nom = cat_vals.get(cont)
                if u_nom is None: continue
                idx_nom = np.searchsorted(umbrales, u_nom)
                if 0 < idx_nom < len(sr_s):
                    col_nom = "#FF7E00" if cat_key=="mala" else "#CC0000"
                    ax_rob.scatter(sr_s[idx_nom], pod_s[idx_nom],
                                   c=col_nom, s=120, marker="D", zorder=7,
                                   edgecolors="white", lw=0.9)
                    ax_rob.annotate(
                        f"NOM {cat_key.replace('_',' ')}",
                        xy=(sr_s[idx_nom], pod_s[idx_nom]),
                        xytext=(sr_s[idx_nom]+0.03, pod_s[idx_nom]+0.025),
                        fontsize=6.5, color=col_nom,
                        path_effects=[pe.withStroke(linewidth=1.5,
                                                    foreground="white")],
                        zorder=8,
                    )

    # Barra de color para la trayectoria
    sm = plt.cm.ScalarMappable(cmap=cmap_tray, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_rob, fraction=0.03, pad=0.02)
    cbar.set_label(f"Umbral normalizado  [{unidad}]", fontsize=7.5)
    import matplotlib.ticker as mticker
    cbar.ax.yaxis.set_major_locator(mticker.FixedLocator([0, 0.2, 0.4, 0.6, 0.8, 1.0]))
    cbar.ax.set_yticklabels(["Bajo","","","","","Alto"], fontsize=7)

    # Leyendas (deduplicadas)
    for ax in [ax_pod, ax_csi, ax_bias]:
        h, l = ax.get_legend_handles_labels()
        seen = {}
        for hi, li in zip(h, l):
            if li not in seen: seen[li] = hi
        ncol = 2 if len(seen) > 5 else 1
        ax.legend(seen.values(), seen.keys(), fontsize=7,
                  loc="best", framealpha=0.88, ncol=ncol)
    ax_rob.legend(fontsize=6.5, loc="lower left",
                  framealpha=0.88, ncol=2)

    fig.suptitle(titulo, fontsize=12, fontweight="bold", y=1.01)
    fig.text(0.01, -0.01,
             "Roebber (2009) Wea. Forecasting 24:749-755  |  "
             "NOM-172-SEMARNAT-2023  |  WRF-Chem vs SINAICA/INECC  |  ICAyCC, UNAM",
             fontsize=6.5, style="italic", color="gray")

    plt.tight_layout()
    os.makedirs(dir_salida, exist_ok=True)
    ruta = os.path.join(dir_salida, nombre_archivo)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK]   {ruta}")


# ──────────────────────────────────────────────────────────────────────────────
# LOGICA DE AGRUPACION
# ──────────────────────────────────────────────────────────────────────────────

def procesar_cont_horizonte(
    datos: pd.DataFrame, cont: str, col_mod: str,
    paso: float, agrupar: str,
    dir_salida: str, suavizar: bool,
) -> pd.DataFrame:

    meta      = META_CONT[cont]
    hor_label = HORIZONTES_LABELS[col_mod]
    df_cont   = datos[datos["contaminante"] == cont].copy()
    if df_cont.empty:
        print(f"[WARN] {cont}: sin datos.")
        return pd.DataFrame()

    df_cont["mes"]       = df_cont["fecha"].dt.month
    df_cont["temporada"] = df_cont["mes"].map(mes_a_temporada)

    ciudades = sorted(df_cont["ciudad"].unique(),
                      key=lambda c: CIUDADES_DOMINIO.index(c)
                      if c in CIUDADES_DOMINIO else 99)

    # Rango automatico sobre observaciones globales limpias
    obs_global = []
    for ciudad in ciudades:
        dfc = df_cont[df_cont["ciudad"] == ciudad]
        o = np.asarray(dfc["max_obs"].values, float)
        m = np.asarray(dfc[col_mod].values, float)
        o, _ = limpiar_serie(o, m, cont)
        obs_global.extend(o.tolist())
    obs_global = np.array(obs_global)
    if len(obs_global) < 10:
        print(f"[WARN] {cont}: datos insuficientes tras QC.")
        return pd.DataFrame()

    u_min, u_max = rango_automatico(obs_global, cont, paso)
    print(f"       Rango barrido: [{u_min:.1f}, {u_max:.1f}] {meta['unidad']}  paso={paso}")

    todos_csv = []
    resultados: Dict = {}
    paleta:     Dict = {}

    if agrupar == "ciudad":
        for ciudad in ciudades:
            dfc = df_cont[df_cont["ciudad"] == ciudad]
            obs = np.asarray(dfc["max_obs"].values, float)
            mod = np.asarray(dfc[col_mod].values, float)
            obs, mod = limpiar_serie(obs, mod, cont)
            df_b = barrido_umbrales(obs, mod, u_min, u_max, paso)
            df_b["ciudad"] = ciudad; df_b["horizonte"] = hor_label
            df_b["agrupacion"] = "ciudad"; df_b["contaminante"] = cont
            todos_csv.append(df_b)
            resultados[ciudad] = df_b
            paleta[ciudad]     = COLORES_CIUDAD.get(ciudad, "#333")
            u_opt, csi_opt = calcular_csi_max(df_b)
            print(f"       {ciudad:<12}: CSI_max={csi_opt:.3f} @ {u_opt:.1f} {meta['unidad']}")

        titulo  = (f"Curva de habilidad — {meta['nombre']}  |  {hor_label}\n"
                   f"Agrupacion por ciudad  |  "
                   f"Rango: [{u_min:.0f}, {u_max:.0f}] {meta['unidad']}")
        nombre  = f"curva_habilidad_{cont}_{hor_label.replace(' ','')}_ciudad.png"

    elif agrupar == "temporada":
        for temp, cfg in TEMPORADAS.items():
            dft = df_cont[df_cont["temporada"] == temp]
            if dft.empty: continue
            obs = np.asarray(dft["max_obs"].values, float)
            mod = np.asarray(dft[col_mod].values, float)
            obs, mod = limpiar_serie(obs, mod, cont)
            df_b = barrido_umbrales(obs, mod, u_min, u_max, paso)
            df_b["temporada"] = temp; df_b["horizonte"] = hor_label
            df_b["agrupacion"] = "temporada"; df_b["contaminante"] = cont
            todos_csv.append(df_b)
            resultados[temp] = df_b
            paleta[temp]     = cfg["color"]
            u_opt, csi_opt = calcular_csi_max(df_b)
            print(f"       {temp:<22}: CSI_max={csi_opt:.3f} @ {u_opt:.1f} {meta['unidad']}")

        titulo  = (f"Curva de habilidad — {meta['nombre']}  |  {hor_label}\n"
                   f"Agrupacion por temporada climatica  |  "
                   f"Rango: [{u_min:.0f}, {u_max:.0f}] {meta['unidad']}")
        nombre  = f"curva_habilidad_{cont}_{hor_label.replace(' ','')}_temporada.png"

    else:   # horizonte: todos en un diagrama, una curva por horizonte
        for c_mod, c_lbl in HORIZONTES_LABELS.items():
            obs_all, mod_all = [], []
            for ciudad in ciudades:
                dfc = df_cont[df_cont["ciudad"] == ciudad]
                o = np.asarray(dfc["max_obs"].values, float)
                m = np.asarray(dfc[c_mod].values, float)
                o, m = limpiar_serie(o, m, cont)
                obs_all.extend(o); mod_all.extend(m)
            if len(obs_all) < 10: continue
            obs_a = np.array(obs_all); mod_a = np.array(mod_all)
            df_b = barrido_umbrales(obs_a, mod_a, u_min, u_max, paso)
            df_b["horizonte"] = c_lbl
            df_b["agrupacion"] = "horizonte"; df_b["contaminante"] = cont
            todos_csv.append(df_b)
            resultados[c_lbl] = df_b
            paleta[c_lbl]     = COLORES_HOR.get(c_mod, "#333")
            u_opt, csi_opt = calcular_csi_max(df_b)
            print(f"       {c_lbl}: CSI_max={csi_opt:.3f} @ {u_opt:.1f} {meta['unidad']}")

        titulo  = (f"Curva de habilidad — {meta['nombre']}  |  todos los horizontes\n"
                   f"Agrupacion por horizonte  |  "
                   f"Rango: [{u_min:.0f}, {u_max:.0f}] {meta['unidad']}")
        nombre  = f"curva_habilidad_{cont}_todos_horizonte.png"

    if resultados:
        graficar_curva_habilidad(resultados, cont, titulo,
                                 dir_salida, nombre,
                                 suavizar=suavizar, paleta=paleta)

    return pd.concat(todos_csv, ignore_index=True) if todos_csv else pd.DataFrame()


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="curva_habilidad.py",
        description=(
            "Analiza la sensibilidad del pronostico WRF-Chem a umbrales "
            "alternativos de concentracion mediante barrido continuo de niveles. "
            "Genera curvas POD/SR/CSI/TSS/BIAS vs umbral y la trayectoria de "
            "habilidad en el espacio de Roebber (2009)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--entrada", "-i", default="combinado/ajustados")
    p.add_argument("--salida",  "-o", default="resultados_curva_habilidad")
    p.add_argument("--cont", "-C", nargs="+", default=None,
                   choices=list(META_CONT.keys()),
                   help="Contaminante(s) a procesar (default: todos)")
    p.add_argument("--horizonte", "-H", default="todos",
                   help="Horizonte: 24h, 48h, 72h o 'todos' (default: todos)")
    p.add_argument("--agrupar", "-g", default="ciudad",
                   choices=["ciudad", "temporada", "horizonte"],
                   help="Agrupacion de las curvas (default: ciudad)")
    p.add_argument("--paso", type=float, default=5.0,
                   help="Paso del barrido en unidades nativas (default: 5.0)")
    p.add_argument("--min-eventos", type=int, default=MIN_EVENTOS,
                   help=f"Minimo de eventos obs por punto (default: {MIN_EVENTOS})")
    p.add_argument("--ciudades", "-c", nargs="+", default=None)
    p.add_argument("--inicio", default=None, metavar="YYYY-MM")
    p.add_argument("--fin",    default=None, metavar="YYYY-MM")
    p.add_argument("--umbral-pm25", type=float, default=UMBRAL_PM25)
    p.add_argument("--iqr-factor",  type=float, default=IQR_FACTOR)
    p.add_argument("--sin-suavizar", action="store_true")
    p.add_argument("--dpi", type=int, default=DPI)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    global IQR_FACTOR, UMBRAL_PM25, MIN_EVENTOS, DPI
    IQR_FACTOR  = args.iqr_factor
    UMBRAL_PM25 = args.umbral_pm25
    MIN_EVENTOS = args.min_eventos
    DPI         = args.dpi

    # Horizontes a procesar
    if args.horizonte.lower() in ("todos", "all"):
        if args.agrupar == "horizonte":
            cols_mod = ["mod_dia1"]   # el modo horizonte maneja los tres internamente
        else:
            cols_mod = list(HORIZONTES_LABELS.keys())
    else:
        col = HORIZONTES_MAP.get(args.horizonte.lower().replace("+",""))
        if col is None:
            sys.exit(f"[ERROR] Horizonte no reconocido: '{args.horizonte}'")
        cols_mod = [col]

    conts           = args.cont or list(META_CONT.keys())
    ciudades_filtro = parsear_lista_ciudades(args.ciudades)

    print("=" * 66)
    print("  Curva de Habilidad — WRF-Chem / SINAICA")
    print("=" * 66)
    print(f"  Entrada      : {args.entrada}")
    print(f"  Salida       : {args.salida}")
    print(f"  Contaminantes: {conts}")
    print(f"  Horizontes   : {[HORIZONTES_LABELS[c] for c in cols_mod]}")
    print(f"  Agrupacion   : {args.agrupar}")
    print(f"  Paso barrido : {args.paso}")
    print(f"  Ciudades     : {ciudades_filtro or 'todas'}")
    print(f"  QC PM2.5     : techo={UMBRAL_PM25} ug/m3  IQR*k={IQR_FACTOR}")
    print("=" * 66)

    datos = leer_datos(args.entrada, ciudades_filtro,
                       args.inicio, args.fin, None)

    csv_acumulado = []
    for cont in conts:
        if cont not in datos["contaminante"].unique():
            print(f"[WARN] {cont}: sin datos — omitido.")
            continue
        for col_mod in cols_mod:
            hor_label = HORIZONTES_LABELS[col_mod]
            print(f"\n-- {cont}  {hor_label}  agrupacion={args.agrupar} --")
            df_b = procesar_cont_horizonte(
                datos, cont, col_mod, args.paso,
                args.agrupar, args.salida, not args.sin_suavizar,
            )
            if not df_b.empty:
                csv_acumulado.append(df_b)

    if csv_acumulado:
        df_csv = pd.concat(csv_acumulado, ignore_index=True)
        cols_ord = ["contaminante","agrupacion","horizonte",
                    "ciudad","temporada","umbral",
                    "N","N_ev_obs","N_ev_mod","frec_obs","frec_mod",
                    "H","M","F","C","POD","SR","FAR","CSI","TSS","BIAS"]
        cols_pres = [c for c in cols_ord if c in df_csv.columns]
        sort_by   = [c for c in ["contaminante","agrupacion","horizonte",
                                  "ciudad","temporada","umbral"]
                     if c in cols_pres]
        df_csv = df_csv[cols_pres].sort_values(sort_by).round(4)
        os.makedirs(args.salida, exist_ok=True)
        ruta_csv = os.path.join(args.salida, "tabla_barrido_umbrales.csv")
        df_csv.to_csv(ruta_csv, index=False)
        print(f"\n[OK]  CSV: {ruta_csv}  ({len(df_csv)} filas)")

    print(f"\n[DONE] Resultados en '{args.salida}/'")


if __name__ == "__main__":
    main()
