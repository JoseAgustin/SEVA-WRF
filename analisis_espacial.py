#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analisis_espacial.py
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
"""