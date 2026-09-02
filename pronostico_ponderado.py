#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pronostico_ponderado.py
=======================
Combina los tres horizontes de pronóstico WRF-Chem (+24h, +48h, +72h)
mediante pesos derivados del desempeño histórico (RMSE inverso) para
producir un pronóstico ponderado para HOY y MAÑANA, separado por
temporada climática (secas frías, secas calientes, lluvias).

Esquema de horizontes
---------------------
  Archivo eval_<CONT>_<Ciudad>_FECHA.csv:
    mod_dia1 → run de (FECHA-1) → +24h → válido para FECHA
    mod_dia2 → run de (FECHA-2) → +48h → válido para FECHA
    mod_dia3 → run de (FECHA-3) → +72h → válido para FECHA

  HOY  (T)   : P = w1·dia1(T) + w2·dia2(T) + w3·dia3(T)
  MAÑANA(T+1): P = w1'·dia1(T+1) + w2'·dia2(T)   [pesos renormalizados]

Pesos
-----
  rmse_inv  : w_i = (1/RMSE_i) / Σ(1/RMSE_j)  [default]
  igual     : w_i = 1/3
  manual    : --w1 --w2 --w3  (se renormalizan)

Temporadas (centro de México)
------------------------------
  Secas frías      : nov–feb
  Secas calientes  : mar–may  (máximo O₃)
  Lluvias          : jun–oct

Salidas
-------
  pesos_historicos.csv
  pronostico_ponderado_HOY.csv
  pronostico_ponderado_MANANA.csv
  serie_ponderada_historica.csv
  pesos_por_temporada.png
  pronostico_ponderado_<CONT>_serie.png
  pronostico_ponderado_scatter.png

Uso
---
  python3 pronostico_ponderado.py
  python3 pronostico_ponderado.py --pesos rmse_inv --fecha 2026-06-30
  python3 pronostico_ponderado.py --pesos manual --w1 0.5 --w2 0.3 --w3 0.2
  python3 pronostico_ponderado.py --sin-temporada --cont O3

Autor  : Pipeline ddsinaica / WRF-Chem — ICAyCC, UNAM
Versión: 1.0.0 (2026-08)
"""

import argparse, glob, os, re, sys
from datetime import date, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

# ── Configuración global ──────────────────────────────────────────────────────

META_CONT = {
    "O3":   {"nombre": "Ozono (O3)",   "unidad": "ppbv",  "color": "#1f77b4"},
    "PM10": {"nombre": "PM10",          "unidad": "ug/m3", "color": "#d62728"},
    "PM25": {"nombre": "PM2.5",         "unidad": "ug/m3", "color": "#ff7f0e"},
    "SO2":  {"nombre": "SO2",           "unidad": "ppbv",  "color": "#9467bd"},
}

NOM172 = {
    "mala":     {"O3":135.0,"PM10":132.0,"PM25": 79.0,"SO2":185.0},
    "muy_mala": {"O3":175.0,"PM10":213.0,"PM25":130.0,"SO2":304.0},
}

CIUDADES_DOMINIO = [
    "CDMX","Cuernavaca","Pachuca","Puebla",
    "SJdelRio","Tlaxcala","Toluca","Tula",
]
_CIUDADES_NORM = {c.lower(): c for c in CIUDADES_DOMINIO}
_CIUDADES_NORM.update({"sjdelrio":"SJdelRio","cdmx":"CDMX"})

TEMPORADAS = {
    "Secas frias":     {"meses":{11,12,1,2}, "color":"#4A90D9"},
    "Secas calientes": {"meses":{3,4,5},      "color":"#F5A623"},
    "Lluvias":         {"meses":{6,7,8,9,10}, "color":"#2d8a2d"},
}

COLS_MOD    = ["mod_dia1","mod_dia2","mod_dia3"]
HOR_LABELS  = {c: f"+{24*(i+1)} h" for i,c in enumerate(COLS_MOD)}
LIMITES     = {
    "O3":  (0,300,0,300),"PM10":(0,1000,0,1000),
    "PM25":(0,500,0,500),"SO2": (0,1000,0,1000),
}
UMBRAL_PM25 = 500.0
IQR_FACTOR  = 3.0
MIN_DIAS    = 10
DPI         = 150

# ── Utilidades ────────────────────────────────────────────────────────────────

def norm_ciudad(n):
    return _CIUDADES_NORM.get(n.strip().lower())

def parsear_ciudades(valor):
    if not valor: return None
    crudos = []
    for item in ([valor] if isinstance(valor,str) else valor):
        crudos.extend(t for t in re.split(r"[,\s]+",item.strip()) if t)
    can, inv = [], []
    for c in crudos:
        cn = norm_ciudad(c)
        if cn is None: inv.append(c)
        elif cn not in can: can.append(cn)
    if inv: sys.exit(f"[ERROR] Ciudades no reconocidas: {inv}")
    return can

def mes_temp(mes):
    for n,cfg in TEMPORADAS.items():
        if mes in cfg["meses"]: return n
    return "Lluvias"

def parsear_nombre(ruta):
    m = re.match(r"^eval_([A-Za-z0-9]+)_(.+)_(\d{4}-\d{2}-\d{2})$",
                 Path(ruta).stem)
    if not m: return None
    return m.group(1).upper(), norm_ciudad(m.group(2)) or m.group(2), m.group(3)

# ── Control de calidad ────────────────────────────────────────────────────────

def limpiar(df, cont):
    lim = LIMITES.get(cont,(0,1e9,0,1e9))
    mo,xo,mm,xm = lim
    if cont=="PM25": xo=UMBRAL_PM25
    df = df[df["max_obs"].between(mo,xo)].copy()
    for c in COLS_MOD:
        df = df[df[c].between(mm,xm)]
    if cont=="PM25" and len(df)>=4:
        q1,q3 = df["max_obs"].quantile([0.25,0.75])
        iqr = q3-q1
        if iqr>0:
            df = df[df["max_obs"].between(q1-IQR_FACTOR*iqr, q3+IQR_FACTOR*iqr)]
    return df

# ── Lectura de datos ──────────────────────────────────────────────────────────

def leer(directorio, ciudades_f, cont_f):
    archivos = sorted(glob.glob(os.path.join(directorio,"eval_*.csv")))
    if not archivos: sys.exit(f"[ERROR] Sin eval_*.csv en '{directorio}'")
    partes = []
    for ruta in archivos:
        meta = parsear_nombre(ruta)
        if not meta: continue
        cont, ciudad, fecha = meta
        if cont not in META_CONT: continue
        if cont_f and cont != cont_f: continue
        if ciudades_f and ciudad not in ciudades_f: continue
        try:
            df = pd.read_csv(ruta, parse_dates=["Fecha"])
        except Exception: continue
        if not {"Fecha","max_obs","mod_dia1","mod_dia2","mod_dia3"}.issubset(df.columns): continue
        df["ciudad"] = ciudad; df["contaminante"] = cont
        partes.append(df[["Fecha","ciudad","contaminante",
                           "max_obs","mod_dia1","mod_dia2","mod_dia3"]])
    if not partes: sys.exit("[ERROR] Ningún archivo coincidió.")
    hist = pd.concat(partes,ignore_index=True).rename(columns={"Fecha":"fecha"})
    hist["fecha"]     = pd.to_datetime(hist["fecha"])
    hist["mes"]       = hist["fecha"].dt.month
    hist["temporada"] = hist["mes"].map(mes_temp)
    # QC por contaminante
    partes2 = []
    for cont in hist["contaminante"].unique():
        partes2.append(limpiar(hist[hist["contaminante"]==cont], cont))
    hist = pd.concat(partes2,ignore_index=True)
    print(f"[INFO] Histórico: {len(hist):,} registros  "
          f"({hist['fecha'].dt.date.nunique()} días, "
          f"{hist['ciudad'].nunique()} ciudades, "
          f"{hist['contaminante'].nunique()} contaminantes)")
    print(f"[INFO] Período  : {hist['fecha'].min().date()} → {hist['fecha'].max().date()}")
    return hist.sort_values(["fecha","ciudad","contaminante"]).reset_index(drop=True)

# ── Cálculo de pesos ──────────────────────────────────────────────────────────

def metricas_hor(obs, mod):
    mask = np.isfinite(obs) & np.isfinite(mod)
    o, m = obs[mask], mod[mask]
    n = len(o)
    if n < 3: return {"RMSE":np.nan,"MAE":np.nan,"R":np.nan,"n":n}
    return {"RMSE": float(np.sqrt(np.mean((m-o)**2))),
            "MAE":  float(np.mean(np.abs(m-o))),
            "R":    float(np.corrcoef(o,m)[0,1]) if np.std(o)>1e-10 else np.nan,
            "n": n}

def calcular_pesos(hist, modo, w_manual, por_temporada):
    filas = []
    for cont in sorted(hist["contaminante"].unique()):
        grupos = ([(cont,t) for t in TEMPORADAS] if por_temporada
                  else [(cont,"Global")])
        for cont_, temp in grupos:
            df_g = hist[hist["contaminante"]==cont_]
            if temp != "Global": df_g = df_g[df_g["temporada"]==temp]
            if len(df_g) < MIN_DIAS: continue
            obs = df_g["max_obs"].values.astype(float)
            met = {c: metricas_hor(obs, df_g[c].values.astype(float))
                   for c in COLS_MOD}
            fila = {"contaminante": cont_, "temporada": temp}
            for i,c in enumerate(COLS_MOD,1):
                for k,v in met[c].items():
                    fila[f"{k}_d{i}"] = v
            # Pesos
            if modo == "igual":
                fila["w1"]=fila["w2"]=fila["w3"]=1/3
            elif modo == "manual" and w_manual:
                t = sum(w_manual)
                fila["w1"],fila["w2"],fila["w3"] = [w/t for w in w_manual]
            else:   # rmse_inv
                rmses = [fila.get(f"RMSE_d{i}", np.nan) for i in range(1,4)]
                if any(np.isnan(r) or r<=0 for r in rmses):
                    fila["w1"]=fila["w2"]=fila["w3"]=1/3
                else:
                    inv=[1/r for r in rmses]; t=sum(inv)
                    fila["w1"],fila["w2"],fila["w3"]=[v/t for v in inv]
            filas.append(fila)
    return pd.DataFrame(filas)

def get_pesos(df_p, cont, temp, por_temporada):
    key = temp if por_temporada else "Global"
    row = df_p[(df_p["contaminante"]==cont)&(df_p["temporada"]==key)]
    if row.empty: return 1/3, 1/3, 1/3
    r = row.iloc[0]
    return float(r["w1"]), float(r["w2"]), float(r["w3"])

def incertidumbre(hist, df_p, cont, temp, w1, w2, w3, por_temporada):
    df_g = hist[hist["contaminante"]==cont]
    if por_temporada: df_g = df_g[df_g["temporada"]==temp]
    if len(df_g)<MIN_DIAS: return np.nan
    obs = df_g["max_obs"].values.astype(float)
    pond = (w1*df_g["mod_dia1"].values.astype(float) +
            w2*df_g["mod_dia2"].values.astype(float) +
            w3*df_g["mod_dia3"].values.astype(float))
    err = pond - obs
    return float(np.nanstd(err))

# ── Pronóstico HOY y MAÑANA ───────────────────────────────────────────────────

def generar(hist, df_p, fecha_hoy, por_temporada):
    ts_hoy = pd.Timestamp(fecha_hoy)
    ts_man = pd.Timestamp(fecha_hoy + timedelta(days=1))
    fh, fm = [], []
    ciudades = sorted(hist["ciudad"].unique(),
                      key=lambda c: CIUDADES_DOMINIO.index(c)
                      if c in CIUDADES_DOMINIO else 99)
    for cont in sorted(hist["contaminante"].unique()):
        meta = META_CONT[cont]
        for ciudad in ciudades:
            df_cc = hist[(hist["contaminante"]==cont)&(hist["ciudad"]==ciudad)]

            # HOY
            rh = df_cc[df_cc["fecha"]==ts_hoy]
            if not rh.empty:
                r   = rh.iloc[0]
                tmp = mes_temp(fecha_hoy.month)
                w1,w2,w3 = get_pesos(df_p,cont,tmp,por_temporada)
                d1,d2,d3 = float(r["mod_dia1"]),float(r["mod_dia2"]),float(r["mod_dia3"])
                disp = [np.isfinite(v) for v in [d1,d2,d3]]
                ws   = [w1,w2,w3]
                tot  = sum(w for w,ok in zip(ws,disp) if ok)
                wadj = [w/tot if ok else 0 for w,ok in zip(ws,disp)]
                pond = sum(w*v for w,v,ok in zip(wadj,[d1,d2,d3],disp) if ok)
                sig  = incertidumbre(hist,df_p,cont,tmp,w1,w2,w3,por_temporada)
                fh.append({"fecha":fecha_hoy.isoformat(),"ciudad":ciudad,
                            "contaminante":cont,"unidad":meta["unidad"],
                            "temporada":tmp,
                            "mod_dia1":round(d1,3),"mod_dia2":round(d2,3),
                            "mod_dia3":round(d3,3),
                            "w1":round(wadj[0],4),"w2":round(wadj[1],4),
                            "w3":round(wadj[2],4),
                            "horizontes_usados":sum(disp),
                            "mod_ponderado_hoy":round(pond,3),
                            "incertidumbre_1s":round(sig,3) if np.isfinite(sig) else None,
                            "obs_hoy":round(float(r["max_obs"]),3)
                                      if np.isfinite(float(r["max_obs"])) else None})

            # MAÑANA
            rm  = df_cc[df_cc["fecha"]==ts_man]
            rh2 = df_cc[df_cc["fecha"]==ts_hoy]
            tmp_m = mes_temp(ts_man.month)
            w1,w2,_ = get_pesos(df_p,cont,tmp_m,por_temporada)
            d1m = float(rm.iloc[0]["mod_dia1"]) if not rm.empty and np.isfinite(float(rm.iloc[0]["mod_dia1"])) else np.nan
            d2m = float(rh2.iloc[0]["mod_dia2"]) if not rh2.empty and np.isfinite(float(rh2.iloc[0]["mod_dia2"])) else np.nan
            disp_m = [np.isfinite(d1m), np.isfinite(d2m)]
            ws_m   = [w1,w2]
            tot_m  = sum(w for w,ok in zip(ws_m,disp_m) if ok)
            wadj_m = [w/tot_m if ok else 0 for w,ok in zip(ws_m,disp_m)] if tot_m>0 else [0.5,0.5]
            pond_m = sum(w*v for w,v,ok in zip(wadj_m,[d1m,d2m],disp_m) if ok)
            sig_m  = incertidumbre(hist,df_p,cont,tmp_m,w1,w2,0,por_temporada)
            nota   = ("dia1+dia2" if np.isfinite(d1m) else
                      "solo dia2 (WRF hoy pendiente)" if np.isfinite(d2m) else
                      "sin datos")
            fm.append({"fecha":(fecha_hoy+timedelta(1)).isoformat(),
                       "ciudad":ciudad,"contaminante":cont,
                       "unidad":meta["unidad"],"temporada":tmp_m,
                       "mod_dia1_man":round(d1m,3) if np.isfinite(d1m) else None,
                       "mod_dia2_man":round(d2m,3) if np.isfinite(d2m) else None,
                       "w1":round(wadj_m[0],4),"w2":round(wadj_m[1],4),
                       "horizontes_usados":sum(disp_m),
                       "mod_ponderado_manana":round(pond_m,3) if np.isfinite(pond_m) else None,
                       "incertidumbre_1s":round(sig_m,3) if np.isfinite(sig_m) else None,
                       "nota":nota})
    return pd.DataFrame(fh), pd.DataFrame(fm)

# ── Serie histórica ponderada ─────────────────────────────────────────────────

def serie_historica(hist, df_p, por_temporada):
    filas = []
    for cont in sorted(hist["contaminante"].unique()):
        for ciudad in sorted(hist["ciudad"].unique()):
            df_cc = hist[(hist["contaminante"]==cont)&(hist["ciudad"]==ciudad)]
            for _,row in df_cc.iterrows():
                tmp = row["temporada"]
                w1,w2,w3 = get_pesos(df_p,cont,tmp,por_temporada)
                d1=float(row["mod_dia1"]); d2=float(row["mod_dia2"]); d3=float(row["mod_dia3"])
                filas.append({"fecha":row["fecha"],"ciudad":ciudad,
                               "contaminante":cont,"temporada":tmp,
                               "obs":float(row["max_obs"]),
                               "pond":w1*d1+w2*d2+w3*d3,
                               "d1":d1,"d2":d2,"d3":d3,
                               "w1":w1,"w2":w2,"w3":w3})
    return pd.DataFrame(filas).sort_values(["contaminante","ciudad","fecha"])

# ── Gráficas ──────────────────────────────────────────────────────────────────

def graficar_serie(df_s, cont, dir_salida):
    meta = META_CONT[cont]
    df_c = df_s[df_s["contaminante"]==cont]
    ciudades = sorted(df_c["ciudad"].unique(),
                      key=lambda c: CIUDADES_DOMINIO.index(c)
                      if c in CIUDADES_DOMINIO else 99)
    n_c = len(ciudades)
    if n_c == 0: return
    fig, axes = plt.subplots(n_c,1,figsize=(14,2.8*n_c+1.5),sharex=True)
    if n_c==1: axes=[axes]
    for ax, ciudad in zip(axes, ciudades):
        df_ci = df_c[df_c["ciudad"]==ciudad].sort_values("fecha")
        fechas = df_ci["fecha"].values
        obs    = df_ci["obs"].values.astype(float)
        pond   = df_ci["pond"].values.astype(float)
        err    = pond - obs; sigma = np.nanstd(err)
        # Fondo temporadas
        prev_t=None; si=0
        for j,t in enumerate(df_ci["temporada"].values):
            if t!=prev_t:
                if prev_t and j>si:
                    ax.axvspan(fechas[si],fechas[j],
                               color=TEMPORADAS[prev_t]["color"],alpha=0.12,zorder=0)
                prev_t=t; si=j
        if prev_t: ax.axvspan(fechas[si],fechas[-1],
                              color=TEMPORADAS[prev_t]["color"],alpha=0.12,zorder=0)
        ax.fill_between(fechas,pond-sigma,pond+sigma,color="#aaa",alpha=0.20,
                        label="±1σ",zorder=1)
        ax.plot(fechas,obs,"k-",lw=1.3,alpha=0.85,label="Observado",zorder=4)
        ax.plot(fechas,pond,color=meta["color"],lw=1.8,label="Ponderado",zorder=5)
        for col,lbl,col_h in [("d1","+24h","#a8d8ea"),("d2","+48h","#f9ca74"),("d3","+72h","#e77f67")]:
            ax.plot(fechas,df_ci[col].values.astype(float),lw=0.7,alpha=0.5,ls=":",label=lbl)
        for cat,cat_v in NOM172.items():
            u=cat_v.get(cont)
            if u is None: continue
            col_n="#FF7E00" if cat=="mala" else "#CC0000"
            ax.axhline(u,color=col_n,lw=1.0,ls="--",alpha=0.75)
            ax.text(fechas[-1],u,f" NOM {cat} ({u})",fontsize=6,color=col_n,va="bottom")
        mask=np.isfinite(obs)&np.isfinite(pond)
        if mask.sum()>2:
            rmse=np.sqrt(np.mean((pond[mask]-obs[mask])**2))
            r=np.corrcoef(obs[mask],pond[mask])[0,1]
            ax.text(0.01,0.97,f"RMSE={rmse:.2f} {meta['unidad']}  R={r:.3f}  σ={sigma:.2f}",
                    transform=ax.transAxes,fontsize=7.5,va="top",ha="left",
                    bbox=dict(fc="white",ec="gray",alpha=0.8,pad=2))
        ax.set_ylabel(f"{ciudad}\n[{meta['unidad']}]",fontsize=8.5)
        ax.grid(True,color="#DDD",lw=0.5); ax.tick_params(labelsize=8)
        if ax==axes[0]: ax.legend(ncol=4,fontsize=7,loc="upper right",framealpha=0.90)
    # Leyenda temporadas
    parches=[mpatches.Patch(color=cfg["color"],alpha=0.6,label=n)
             for n,cfg in TEMPORADAS.items()]
    axes[-1].legend(handles=parches,ncol=3,fontsize=7,loc="upper left",
                    bbox_to_anchor=(0,-0.22),frameon=True)
    axes[-1].set_xlabel("Fecha",fontsize=9)
    fig.suptitle(f"Pronóstico ponderado — {meta['nombre']}\n"
                 "Negro=obs · Color=ponderado · Punteado=horizontes · Banda=±1σ",
                 fontsize=11,fontweight="bold")
    fig.text(0.01,0.005,"WRF-Chem vs SINAICA/INECC | NOM-172-SEMARNAT-2023 | ICAyCC, UNAM",
             fontsize=6.5,style="italic",color="gray")
    plt.tight_layout(rect=[0,0.05,1,0.97])
    os.makedirs(dir_salida,exist_ok=True)
    ruta=os.path.join(dir_salida,f"pronostico_ponderado_{cont}_serie.png")
    fig.savefig(ruta,dpi=DPI,bbox_inches="tight"); plt.close(fig)
    print(f"[OK]   {ruta}")

def graficar_pesos(df_p, dir_salida):
    conts=sorted(df_p["contaminante"].unique())
    temps=sorted(df_p["temporada"].unique())
    fig,axes=plt.subplots(1,len(conts),figsize=(4*len(conts)+1,4.5),sharey=True)
    if len(conts)==1: axes=[axes]
    colores=[  "#1f77b4","#ff7f0e","#d62728"]
    etqs=["+24h (w1)","+48h (w2)","+72h (w3)"]
    for ax,cont in zip(axes,conts):
        meta=META_CONT[cont]; df_c=df_p[df_p["contaminante"]==cont]
        x=np.arange(len(temps)); bottom=np.zeros(len(temps))
        for i,(cw,col,lbl) in enumerate(zip(["w1","w2","w3"],colores,etqs)):
            vals=[]
            for t in temps:
                r=df_c[df_c["temporada"]==t]
                vals.append(float(r[cw].values[0]) if not r.empty else 1/3)
            ax.bar(x,vals,bottom=bottom,color=col,label=lbl,
                   width=0.6,edgecolor="white",lw=0.8)
            for j,(v,b) in enumerate(zip(vals,bottom)):
                ax.text(x[j],b+v/2,f"{v:.2f}",ha="center",va="center",
                        fontsize=8,color="white",fontweight="bold")
            bottom=bottom+np.array(vals)
        ax2=ax.twinx()
        for i,cr in enumerate(["RMSE_d1","RMSE_d2","RMSE_d3"]):
            if cr not in df_c.columns: continue
            rv=[]
            for t in temps:
                r=df_c[df_c["temporada"]==t]
                rv.append(float(r[cr].values[0]) if not r.empty else np.nan)
            ax2.plot(x,rv,color=colores[i],lw=1.8,ls="--",marker="D",ms=5)
        ax2.set_ylabel(f"RMSE [{meta['unidad']}]",fontsize=8)
        ax2.tick_params(labelsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels([t.replace(" ","\n") for t in temps],fontsize=8)
        ax.set_ylim(0,1.02); ax.set_title(meta["nombre"],fontsize=10,fontweight="bold")
        ax.set_xlabel("Temporada",fontsize=9)
        if ax==axes[0]: ax.set_ylabel("Peso del horizonte",fontsize=9)
        if ax==axes[0]: ax.legend(loc="lower left",fontsize=7,framealpha=0.9)
    fig.suptitle("Pesos de ponderación por horizonte y temporada\n"
                 "Barras=pesos · Punteado=RMSE",fontsize=11,fontweight="bold")
    plt.tight_layout()
    os.makedirs(dir_salida,exist_ok=True)
    ruta=os.path.join(dir_salida,"pesos_por_temporada.png")
    fig.savefig(ruta,dpi=DPI,bbox_inches="tight"); plt.close(fig)
    print(f"[OK]   {ruta}")

def graficar_scatter(df_s, dir_salida):
    conts=sorted(df_s["contaminante"].unique())
    fig,axes=plt.subplots(1,len(conts),figsize=(4.5*len(conts),4.5))
    if len(conts)==1: axes=[axes]
    for ax,cont in zip(axes,conts):
        meta=META_CONT[cont]
        df_c=df_s[df_s["contaminante"]==cont].dropna(subset=["obs","pond"])
        obs=df_c["obs"].values.astype(float); pond=df_c["pond"].values.astype(float)
        mask=np.isfinite(obs)&np.isfinite(pond); obs,pond=obs[mask],pond[mask]
        if len(obs)<3: continue
        temps=df_c[mask]["temporada"].values
        for t,cfg in TEMPORADAS.items():
            tm=temps==t
            ax.scatter(obs[tm],pond[tm],c=cfg["color"],s=18,alpha=0.65,label=t,zorder=3)
        vmin=min(obs.min(),pond.min()); vmax=max(obs.max(),pond.max())
        ax.plot([vmin,vmax],[vmin,vmax],"k--",lw=1.3,alpha=0.6,label="1:1")
        sl,inter,r,*_=stats.linregress(obs,pond)
        x_r=np.array([vmin,vmax])
        ax.plot(x_r,sl*x_r+inter,color=meta["color"],lw=1.5,label=f"Regresión (R={r:.3f})")
        rmse=np.sqrt(np.mean((pond-obs)**2))
        ax.text(0.05,0.95,f"R={r:.3f}\nRMSE={rmse:.2f} {meta['unidad']}\nn={len(obs)}",
                transform=ax.transAxes,fontsize=8,va="top",
                bbox=dict(fc="white",ec="gray",alpha=0.85,pad=2))
        ax.set_xlabel(f"Observado [{meta['unidad']}]",fontsize=9)
        ax.set_ylabel(f"Ponderado [{meta['unidad']}]",fontsize=9)
        ax.set_title(meta["nombre"],fontsize=10,fontweight="bold")
        ax.set_aspect("equal"); ax.grid(True,color="#DDD",lw=0.5)
        ax.legend(fontsize=7,loc="lower right")
    fig.suptitle("Validación histórica — Pronóstico ponderado vs Observaciones",
                 fontsize=11,fontweight="bold")
    fig.text(0.01,0.005,"WRF-Chem vs SINAICA/INECC | NOM-172 | ICAyCC, UNAM",
             fontsize=6.5,style="italic",color="gray")
    plt.tight_layout()
    os.makedirs(dir_salida,exist_ok=True)
    ruta=os.path.join(dir_salida,"pronostico_ponderado_scatter.png")
    fig.savefig(ruta,dpi=DPI,bbox_inches="tight"); plt.close(fig)
    print(f"[OK]   {ruta}")

# ── CLI y main ────────────────────────────────────────────────────────────────

def parse_args():
    p=argparse.ArgumentParser(prog="pronostico_ponderado.py",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--entrada","-i",default="combinado/ajustados")
    p.add_argument("--salida", "-o",default="resultados_ponderado")
    p.add_argument("--fecha",  "-f",default=None,metavar="YYYY-MM-DD")
    p.add_argument("--cont",   "-C",nargs="+",default=None,choices=list(META_CONT.keys()))
    p.add_argument("--ciudades","-c",nargs="+",default=None)
    p.add_argument("--pesos",  "-p",default="rmse_inv",
                   choices=["rmse_inv","igual","manual"])
    p.add_argument("--w1",type=float,default=0.5)
    p.add_argument("--w2",type=float,default=0.3)
    p.add_argument("--w3",type=float,default=0.2)
    p.add_argument("--sin-temporada",action="store_true")
    p.add_argument("--umbral-pm25",type=float,default=UMBRAL_PM25)
    p.add_argument("--iqr-factor", type=float,default=IQR_FACTOR)
    p.add_argument("--dpi",type=int,default=DPI)
    return p.parse_args()

def main():
    args=parse_args()
    global UMBRAL_PM25,IQR_FACTOR,DPI
    UMBRAL_PM25=args.umbral_pm25; IQR_FACTOR=args.iqr_factor; DPI=args.dpi
    por_temp = not args.sin_temporada
    cds_f = parsear_ciudades(args.ciudades)
    cont_f = args.cont[0] if args.cont and len(args.cont)==1 else None
    w_manual = (args.w1,args.w2,args.w3) if args.pesos=="manual" else None

    print("="*66)
    print("  Pronóstico Ponderado — WRF-Chem / SINAICA")
    print("="*66)
    print(f"  Entrada       : {args.entrada}")
    print(f"  Salida        : {args.salida}")
    print(f"  Método pesos  : {args.pesos}")
    print(f"  Por temporada : {por_temp}")
    print(f"  Ciudades      : {cds_f or 'todas'}")
    print(f"  Contaminantes : {args.cont or 'todos'}")
    print("="*66)

    hist = leer(args.entrada, cds_f, cont_f)
    if args.cont and len(args.cont)>1:
        hist = hist[hist["contaminante"].isin(args.cont)]

    fecha_hoy = date.fromisoformat(args.fecha) if args.fecha \
                else hist["fecha"].max().date()
    print(f"\n[INFO] HOY    : {fecha_hoy}")
    print(f"[INFO] MAÑANA : {fecha_hoy+timedelta(1)}")

    print("\n[INFO] Calculando pesos históricos...")
    df_p = calcular_pesos(hist, args.pesos, w_manual, por_temp)
    os.makedirs(args.salida,exist_ok=True)
    df_p.round(4).to_csv(os.path.join(args.salida,"pesos_historicos.csv"),index=False)
    print(f"[OK]   pesos_historicos.csv")

    # Mostrar tabla de pesos
    print(f"\n  {'Cont':<6} {'Temporada':<22} {'w1':>7} {'w2':>7} {'w3':>7}"
          f" {'RMSE_d1':>9} {'RMSE_d2':>9} {'RMSE_d3':>9}")
    print("  "+"-"*78)
    for _,r in df_p.iterrows():
        print(f"  {r['contaminante']:<6} {r['temporada']:<22}"
              f" {r['w1']:>7.3f} {r['w2']:>7.3f} {r['w3']:>7.3f}"
              f" {r.get('RMSE_d1',np.nan):>9.3f}"
              f" {r.get('RMSE_d2',np.nan):>9.3f}"
              f" {r.get('RMSE_d3',np.nan):>9.3f}")

    print(f"\n[INFO] Generando pronóstico {fecha_hoy} / {fecha_hoy+timedelta(1)}...")
    df_hoy, df_man = generar(hist, df_p, fecha_hoy, por_temp)
    df_hoy.to_csv(os.path.join(args.salida,"pronostico_ponderado_HOY.csv"),index=False)
    df_man.to_csv(os.path.join(args.salida,"pronostico_ponderado_MANANA.csv"),index=False)
    print(f"[OK]   pronostico_ponderado_HOY.csv    ({len(df_hoy)} filas)")
    print(f"[OK]   pronostico_ponderado_MANANA.csv ({len(df_man)} filas)")

    if not df_hoy.empty:
        print(f"\n  HOY ({fecha_hoy}):")
        cols=["ciudad","contaminante","mod_ponderado_hoy","incertidumbre_1s",
              "obs_hoy","w1","w2","w3"]
        print(df_hoy[[c for c in cols if c in df_hoy.columns]].to_string(index=False))

    if not df_man.empty:
        print(f"\n  MAÑANA ({fecha_hoy+timedelta(1)}):")
        cols=["ciudad","contaminante","mod_ponderado_manana","incertidumbre_1s",
              "horizontes_usados","nota"]
        print(df_man[[c for c in cols if c in df_man.columns]].to_string(index=False))

    print("\n[INFO] Construyendo serie histórica ponderada...")
    df_s = serie_historica(hist, df_p, por_temp)
    cols_num = df_s.select_dtypes(include=[np.number]).columns
    df_s[cols_num] = df_s[cols_num].round(3)
    df_s.to_csv(os.path.join(args.salida,"serie_ponderada_historica.csv"),index=False)
    print(f"[OK]   serie_ponderada_historica.csv")

    print("\n[INFO] Generando gráficas...")
    for cont in sorted(hist["contaminante"].unique()):
        graficar_serie(df_s, cont, args.salida)
    graficar_pesos(df_p, args.salida)
    graficar_scatter(df_s, args.salida)

    print(f"\n[DONE] Resultados en '{args.salida}/'")

if __name__=="__main__":
    main()
